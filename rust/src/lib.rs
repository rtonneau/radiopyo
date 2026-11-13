//! `radiopyo_rs._core`: Python binding of radkin, the Rust port of radiopyo.
//!
//! The Python side (`radiopyo.rust`) builds a [`Cell`] from a `UnitCell`'s
//! merged `ConfigDict`, radiopyo's Boscolo CSV and its current state
//! (temperature, beam, LET), checks [`Cell::system`] against its own
//! matrices (parity guard), then calls [`Cell::run`]. The returned dict holds
//! what `ODEResult.append_sol` needs, one entry per solved block.

use std::sync::OnceLock;

use numpy::ndarray::Array2;
use numpy::{IntoPyArray, PyArray1};
use pyo3::create_exception;
use pyo3::exceptions::PyException;
use pyo3::prelude::*;
use pyo3::types::{PyBool, PyDict, PyFloat, PyInt, PyList, PyString, PyTuple};
use rayon::prelude::*;
use radkin::RadkinError;
use radkin::config::Config;
use radkin::gvalues::GValuesBoscolo;
use radkin::physics::Beam;
use radkin::result::OdeResult;
use radkin::solver::JacobianMode;
use radkin::unit_cell::{ChunkOptions, RunOptions, UnitCell};

create_exception!(
    radiopyo_rs._core,
    RadkinRustError,
    PyException,
    "radkin error; args are (kind, message), kind naming the RadkinError variant."
);

/// Kind of a radkin error, mapped to a radiopyo exception on the Python side.
fn error_kind(e: &RadkinError) -> &'static str {
    match e {
        RadkinError::Parse { .. } => "parse",
        RadkinError::Config(_) => "config",
        RadkinError::Io { .. } => "io",
        RadkinError::Beam(_) => "beam",
        RadkinError::ShortBeamOn(_) => "short_beam_on",
        RadkinError::NoBeamDefined => "no_beam",
        RadkinError::SolverFailed { .. } => "solver_failed",
        RadkinError::GValue(_) => "gvalue",
        RadkinError::Value(_) => "value",
    }
}

fn to_py_err(e: RadkinError) -> PyErr {
    RadkinRustError::new_err((error_kind(&e), e.to_string()))
}

fn value_err(msg: impl Into<String>) -> PyErr {
    RadkinRustError::new_err(("value", msg.into()))
}

// ------------------------------- Logging -------------------------------- //

/// Forwards radkin's `log` records to Python's `radiopyo.radkin` logger
/// (a child of `radiopyo`, so radiopyo's handlers print them).
struct PyLogger;

static LOGGER: PyLogger = PyLogger;
static PY_LOGGER: OnceLock<Py<PyAny>> = OnceLock::new();

impl log::Log for PyLogger {
    fn enabled(&self, metadata: &log::Metadata) -> bool {
        metadata.level() <= log::max_level()
    }

    fn log(&self, record: &log::Record) {
        if !self.enabled(record.metadata()) {
            return;
        }
        let Some(logger) = PY_LOGGER.get() else {
            return;
        };
        let level: u8 = match record.level() {
            log::Level::Error => 40,
            log::Level::Warn => 30,
            log::Level::Info => 20,
            log::Level::Debug | log::Level::Trace => 10,
        };
        let msg = record.args().to_string();
        Python::attach(|py| {
            // A failing Python logger must never abort a solve.
            let _ = logger.call_method1(py, "log", (level, msg));
        });
    }

    fn flush(&self) {}
}

fn init_logging(py: Python<'_>) -> PyResult<()> {
    let logger = py
        .import("logging")?
        .call_method1("getLogger", ("radiopyo.radkin",))?;
    let _ = PY_LOGGER.set(logger.unbind());
    // Already set (module re-imported in another interpreter): keep it.
    if log::set_logger(&LOGGER).is_ok() {
        log::set_max_level(log::LevelFilter::Info);
    }
    Ok(())
}

// --------------------------- Python -> TOML ---------------------------- //

/// Convert a Python value (dict/list/str/bool/int/float) into TOML.
fn to_toml(obj: &Bound<'_, PyAny>, path: &str) -> PyResult<toml::Value> {
    if obj.is_instance_of::<PyBool>() {
        return Ok(toml::Value::Boolean(obj.extract()?));
    }
    if obj.is_instance_of::<PyInt>() {
        return Ok(match obj.extract::<i64>() {
            Ok(i) => toml::Value::Integer(i),
            Err(_) => toml::Value::Float(obj.extract::<f64>()?),
        });
    }
    if obj.is_instance_of::<PyFloat>() {
        return Ok(toml::Value::Float(obj.extract()?));
    }
    if obj.is_instance_of::<PyString>() {
        return Ok(toml::Value::String(obj.extract()?));
    }
    if let Ok(d) = obj.cast::<PyDict>() {
        return Ok(toml::Value::Table(dict_to_table(d, path)?));
    }
    if obj.is_instance_of::<PyList>() || obj.is_instance_of::<PyTuple>() {
        let mut out = Vec::new();
        for (i, item) in obj.try_iter()?.enumerate() {
            out.push(to_toml(&item?, &format!("{path}[{i}]"))?);
        }
        return Ok(toml::Value::Array(out));
    }
    // numpy scalars and other numbers.
    if let Ok(f) = obj.extract::<f64>() {
        return Ok(toml::Value::Float(f));
    }
    Err(value_err(format!(
        "config value at '{path}' cannot be passed to radkin: {}",
        obj.repr()?
    )))
}

fn dict_to_table(d: &Bound<'_, PyDict>, path: &str) -> PyResult<toml::Table> {
    let mut t = toml::Table::new();
    for (k, v) in d.iter() {
        let key: String = k.extract()?;
        if v.is_none() {
            continue;
        }
        let sub = if path.is_empty() {
            key.clone()
        } else {
            format!("{path}.{key}")
        };
        t.insert(key, to_toml(&v, &sub)?);
    }
    Ok(t)
}

fn get_opt<'py, T: FromPyObjectOwned<'py>>(
    d: &Bound<'py, PyDict>,
    key: &str,
) -> PyResult<Option<T>> {
    match d.get_item(key)? {
        None => Ok(None),
        Some(v) if v.is_none() => Ok(None),
        Some(v) => Ok(Some(v.extract::<T>().map_err(|e| {
            let e: PyErr = e.into();
            value_err(format!("option '{key}': {e}"))
        })?)),
    }
}

// ------------------------------- Options -------------------------------- //

/// Run and chunk options, owned (no Python object): safe to move to a thread.
#[derive(Debug, Clone)]
struct Job {
    mode: Mode,
    t_span: (f64, f64),
    run: RunOptions,
    chunk: ChunkOptions,
    y0: Option<Vec<f64>>,
}

#[derive(Debug, Clone, Copy)]
enum Mode {
    Run,
    RunTillMaxDose,
    Chunked,
    SinglePulse,
}

const OPTION_KEYS: [&str; 17] = [
    "atol",
    "rtol",
    "max_step",
    "first_step",
    "use_beam",
    "use_jac",
    "jacobian",
    "t_eval",
    "y0",
    "max_step_size_on",
    "max_step_size_off",
    "max_step_size_trans",
    "max_step_size_final",
    "trans_time",
    "points_per_chunk_on",
    "points_per_chunk_off",
    "n_threads",
];

fn parse_job(mode: &str, t_span: (f64, f64), o: &Bound<'_, PyDict>) -> PyResult<Job> {
    for k in o.keys() {
        let k: String = k.extract()?;
        if !OPTION_KEYS.contains(&k.as_str()) {
            return Err(value_err(format!("unknown radkin option '{k}'")));
        }
    }
    let mode = match mode {
        "run" => Mode::Run,
        "run_till_max_dose" => Mode::RunTillMaxDose,
        "chunked" => Mode::Chunked,
        "single_pulse" => Mode::SinglePulse,
        m => return Err(value_err(format!("unknown radkin mode '{m}'"))),
    };
    let jacobian = match get_opt::<String>(o, "jacobian")?.as_deref() {
        None | Some("dense") => JacobianMode::Dense,
        Some("sparse") => JacobianMode::Sparse,
        Some(j) => return Err(value_err(format!("unknown jacobian mode '{j}'"))),
    };
    let run = RunOptions {
        atol: get_opt(o, "atol")?,
        rtol: get_opt(o, "rtol")?,
        max_step: get_opt(o, "max_step")?,
        first_step: get_opt(o, "first_step")?,
        use_beam: get_opt(o, "use_beam")?.unwrap_or(true),
        use_jac: get_opt(o, "use_jac")?.unwrap_or(true),
        jacobian,
        t_eval: get_opt(o, "t_eval")?,
    };
    let y0: Option<Vec<f64>> = get_opt(o, "y0")?;
    let chunk = ChunkOptions {
        run: RunOptions {
            t_eval: None,
            max_step: None,
            ..run.clone()
        },
        max_step_size_on: get_opt(o, "max_step_size_on")?,
        max_step_size_off: get_opt(o, "max_step_size_off")?,
        max_step_size_trans: get_opt(o, "max_step_size_trans")?,
        max_step_size_final: get_opt(o, "max_step_size_final")?,
        trans_time: get_opt(o, "trans_time")?,
        points_per_chunk_on: get_opt(o, "points_per_chunk_on")?,
        points_per_chunk_off: get_opt(o, "points_per_chunk_off")?,
        y0: y0.clone(),
    };
    Ok(Job {
        mode,
        t_span,
        run,
        chunk,
        y0,
    })
}

/// Run `job` on a copy of `cell` (the cell itself is never modified).
fn execute(cell: &UnitCell, job: &Job) -> radkin::Result<OdeResult> {
    let mut cell = cell.clone();
    match job.mode {
        Mode::Run => cell.run(job.t_span, &job.run, job.y0.as_deref()),
        Mode::RunTillMaxDose => {
            if job.y0.is_some() {
                return Err(RadkinError::Value(
                    "y0 is not available with run_till_max_dose".into(),
                ));
            }
            cell.run_till_max_dose(job.t_span.0, &job.run)
        }
        Mode::Chunked => cell.chunked_run(job.t_span, &job.chunk),
        Mode::SinglePulse => cell.single_pulse_run(job.t_span, &job.chunk),
    }
}

// ------------------------------- Results -------------------------------- //

fn result_to_py<'py>(py: Python<'py>, res: &OdeResult) -> PyResult<Bound<'py, PyDict>> {
    let out = PyDict::new(py);
    let t = res.t();
    let y = res.y();
    let n_dyn = res.mat.n_dyn;
    let mut flat = vec![0.0; n_dyn * t.len()];
    for (j, row) in y.iter().enumerate() {
        for (i, v) in row.iter().enumerate() {
            flat[i * t.len() + j] = *v;
        }
    }
    let y = Array2::from_shape_vec((n_dyn, t.len()), flat).map_err(|e| value_err(e.to_string()))?;
    out.set_item("t", PyArray1::from_slice(py, t))?;
    out.set_item("y", y.into_pyarray(py))?;
    out.set_item("n_points", res.n_points.clone())?;
    out.set_item("nfev", res.nfev.clone())?;
    out.set_item("njev", res.njev.clone())?;
    out.set_item("nst", res.nst.clone())?;
    out.set_item("exec_time", res.exec_time.clone())?;
    out.set_item("step_dose_rates", res.step_dose_rates.clone())?;
    out.set_item("burst_dose_rates", res.burst_dose_rates.clone())?;
    let spans: Vec<(f64, f64, f64)> = res
        .short_on_spans
        .iter()
        .map(|s| (s.start, s.stop, s.dose_rate))
        .collect();
    out.set_item("short_on_spans", spans)?;
    out.set_item("thinned", res.thinned)?;
    Ok(out)
}

// --------------------------------- Cell --------------------------------- //

/// A radkin `UnitCell` built from a radiopyo configuration and state.
#[pyclass(frozen, module = "radiopyo_rs._core")]
struct Cell {
    inner: UnitCell,
}

#[pymethods]
impl Cell {
    /// Build the cell.
    ///
    /// Args:
    ///     config: merged radiopyo ``ConfigDict`` (includes resolved; its
    ///         ``includes`` and ``beam`` sections are ignored).
    ///     gvalues_csv: text of radiopyo's Boscolo CSV.
    ///     temperature: [degC].
    ///     beam: ``None`` or ``(type, label, params, LET)`` with ``params``
    ///         the beam's ``as_dict()``.
    #[new]
    #[pyo3(signature = (config, gvalues_csv, temperature, beam=None))]
    fn new(
        config: &Bound<'_, PyDict>,
        gvalues_csv: &str,
        temperature: f64,
        beam: Option<(String, String, Bound<'_, PyDict>, Option<f64>)>,
    ) -> PyResult<Self> {
        let mut table = dict_to_table(config, "")?;
        table.remove("includes");
        table.remove("beam");
        let cfg = Config::from_table(table, None).map_err(to_py_err)?;
        let gv =
            GValuesBoscolo::from_csv_str(gvalues_csv, "radiopyo Boscolo CSV").map_err(to_py_err)?;
        let mut cell = UnitCell::from_config_with_gvalues(&cfg, gv).map_err(to_py_err)?;
        cell.set_temperature(temperature).map_err(to_py_err)?;
        if let Some((kind, label, params, let_)) = beam {
            let b = Beam::from_table(&kind, &dict_to_table(&params, "beam")?).map_err(to_py_err)?;
            match let_ {
                Some(l) => cell.set_beam(b, &label, Some(l)).map_err(to_py_err)?,
                None => {
                    return Err(value_err(
                        "the beam has no LET: set one (UnitCell.change_LET) first",
                    ));
                }
            }
        }
        Ok(Self { inner: cell })
    }

    /// The system radkin built, for the parity guard: species labels (dynamic
    /// then fixed, by index), reaction strings (global index order: k pairs
    /// first, then k-reactions, radiolytic, enzymatic), k values, radiolytic
    /// yields [mol/L/Gy], enzymatic k and Michaelis constants, fixed
    /// concentrations [mol/L] and initial values [µmol/L].
    fn system<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyDict>> {
        let c = &self.inner;
        let out = PyDict::new(py);
        let mut species = c.env.species.sorted_dyn_labels();
        species.extend(c.env.species.sorted_cst_labels());
        out.set_item("species", species)?;
        out.set_item("n_dyn", c.mat.n_dyn)?;
        out.set_item("reactions", c.env.reaction_strings())?;
        out.set_item("k_values", c.mat.k.k_values.clone())?;
        out.set_item("g_yields", c.mat.radiolytic.yields.clone())?;
        out.set_item("ez_k", c.mat.enzymatic.k_values.clone())?;
        out.set_item("ez_km", c.mat.enzymatic.k_michaelis.clone())?;
        out.set_item("cst_sp", c.mat.cst_sp.clone())?;
        out.set_item("initial_values", c.env.initial_values())?;
        out.set_item("temperature", c.env.temperature)?;
        Ok(out)
    }

    /// Run the cell (the GIL is released while solving).
    ///
    /// Args:
    ///     mode: ``"run"``, ``"run_till_max_dose"``, ``"chunked"`` or
    ///         ``"single_pulse"``.
    ///     t_span: (start, stop) [s] (``run_till_max_dose`` uses start only).
    ///     options: solver and chunk options (see ``OPTION_KEYS``).
    #[pyo3(signature = (mode, t_span, options=None))]
    fn run<'py>(
        &self,
        py: Python<'py>,
        mode: &str,
        t_span: (f64, f64),
        options: Option<&Bound<'py, PyDict>>,
    ) -> PyResult<Bound<'py, PyDict>> {
        let empty = PyDict::new(py);
        let job = parse_job(mode, t_span, options.unwrap_or(&empty))?;
        let res = py
            .detach(|| execute(&self.inner, &job))
            .map_err(to_py_err)?;
        result_to_py(py, &res)
    }
}

/// One batch item: (cell, mode, t_span, options).
type BatchItem<'py> = (
    Bound<'py, Cell>,
    String,
    (f64, f64),
    Option<Bound<'py, PyDict>>,
);

/// Run several cells in parallel on a rayon pool, GIL released.
///
/// Args:
///     items: list of ``(cell, mode, t_span, options)`` (as ``Cell.run``).
///     n_threads: pool size (``None``: all cores).
///
/// Returns:
///     One entry per item, in input order: the result dict of ``Cell.run``,
///     or the ``RadkinRustError`` instance of a failed item (never raised).
///     Invalid options raise ``ValueError`` before any solve.
#[pyfunction]
#[pyo3(signature = (items, n_threads=None))]
fn run_batch<'py>(
    py: Python<'py>,
    items: Vec<BatchItem<'py>>,
    n_threads: Option<usize>,
) -> PyResult<Vec<Bound<'py, PyAny>>> {
    let empty = PyDict::new(py);
    let mut jobs = Vec::with_capacity(items.len());
    for (_, mode, t_span, opts) in &items {
        jobs.push(parse_job(mode, *t_span, opts.as_ref().unwrap_or(&empty))?);
    }
    let cells: Vec<&UnitCell> = items.iter().map(|(c, ..)| &c.get().inner).collect();
    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(n_threads.unwrap_or(0))
        .build()
        .map_err(|e| value_err(e.to_string()))?;
    let results: Vec<radkin::Result<OdeResult>> = py.detach(|| {
        pool.install(|| {
            cells
                .par_iter()
                .zip(jobs.par_iter())
                .map(|(cell, job)| execute(cell, job))
                .collect()
        })
    });
    results
        .into_iter()
        .map(|r| match r {
            Ok(res) => Ok(result_to_py(py, &res)?.into_any()),
            Err(e) => Ok(to_py_err(e).into_value(py).into_bound(py).into_any()),
        })
        .collect()
}

/// Names of the accepted options (for the Python side's checks).
#[pyfunction]
fn option_keys() -> Vec<&'static str> {
    OPTION_KEYS.to_vec()
}

#[pymodule]
fn _core(m: &Bound<'_, PyModule>) -> PyResult<()> {
    init_logging(m.py())?;
    m.add("__version__", env!("CARGO_PKG_VERSION"))?;
    m.add("RadkinRustError", m.py().get_type::<RadkinRustError>())?;
    m.add_class::<Cell>()?;
    m.add_function(wrap_pyfunction!(option_keys, m)?)?;
    m.add_function(wrap_pyfunction!(run_batch, m)?)?;
    Ok(())
}
