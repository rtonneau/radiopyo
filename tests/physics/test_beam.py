import pytest

from radiopyo.physics import beam


class TestCstTimer():
    timer = beam.Timer.new_constant()  # ==> All float("inf")

    def test_cst(self):
        assert self.timer.duty_cycle() == 1.0


class TestPulsedTimer():
    period = 5
    on_time = 2
    timer = beam.Timer.new_pulsed(5, 2)

    def test_pulsed(self):
        assert self.timer.duty_cycle() == pytest.approx(self.on_time/self.period)

    def test_period_inf(self):
        timer = beam.Timer(period=float("inf"), on_time=4)
        assert timer.duty_cycle() == 0

    def test_on_time_inf(self):
        timer = beam.Timer(period=4, on_time=float("inf"))
        assert timer.duty_cycle() == 1.0


class TestCstBeam():
    dose_rate = 5
    test_beam = beam.ConstantBeam(dose_rate)

    def test_average_dose_rate(self):
        assert self.test_beam.average_dose_rate() == self.dose_rate

    def test_peak_dose_rate(self):
        assert self.test_beam.peak_dose_rate == self.dose_rate

    def test_time_message(self):
        time = 3.0
        assert self.test_beam.at(3).time == pytest.approx(time)
        assert self.test_beam.at(3).dose_rate() == self.dose_rate


class TestPulsedBeam():
    dose_rate = 5
    period = 5
    on_time = 2
    test_beam = beam.PulsedBeam(dose_rate=dose_rate,
                                period=period,
                                on_time=on_time
                                )

    def test_average_dose_rate(self):
        assert self.test_beam.average_dose_rate() == self.dose_rate

    def test_peak_dose_rate(self):
        assert self.test_beam.peak_dose_rate == pytest.approx(
            self.dose_rate*self.period/self.on_time)

    def test_time_message_at_off(self):
        time = 3.0
        assert self.test_beam.at(time).time == pytest.approx(time)
        assert self.test_beam.at(time).dose_rate() == 0.0

    def test_time_message_at_on(self):
        time = 1.0
        assert self.test_beam.at(time).time == pytest.approx(time)
        assert self.test_beam.at(time).dose_rate(
        ) == self.test_beam.peak_dose_rate


class TestPulsedBeamMaxed(TestPulsedBeam):
    dose_rate = 5
    period = 5
    on_time = 2
    max_dose = (dose_rate * period) * 1.5
    test_beam = beam.PulsedBeam(dose_rate=dose_rate,
                                period=period,
                                on_time=on_time,
                                max_dose=max_dose,
                                )

    def test_max_dose_reached(self):
        time = 20
        assert self.test_beam.at(time).time == pytest.approx(time)
        assert self.test_beam.at(time).dose_rate() == 0.0

    def test_max_dose_not_reached(self):
        assert self.test_beam.at(1.5).dose_rate(
        ) == self.test_beam.peak_dose_rate
        assert self.test_beam.at(3).dose_rate() == 0.0
        assert self.test_beam.at(5.5).dose_rate(
        ) == self.test_beam.peak_dose_rate
        assert self.test_beam.at(6.5).dose_rate() == 0.0
        assert self.test_beam.at(11).dose_rate() == 0.0
