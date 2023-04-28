import pytest

from radiopyo.chemistry import radiolytic


def test_conversion():
    assert radiolytic.ge_to_kr(2.80) == pytest.approx(2.9020e-07)
    assert radiolytic.ge_to_kr(0.62) == pytest.approx(6.4258e-08)
    assert radiolytic.ge_to_kr(0.47) == pytest.approx(4.8712e-08)
    assert radiolytic.ge_to_kr(0.73) == pytest.approx(7.5659e-08)
