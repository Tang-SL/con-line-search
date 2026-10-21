import numpy as np
import pytest
from conlines.reliability import (component_properties, detection_reliability, neighbour_counts,
                                  identification_probabilities, chance_match_probability, anchor_statistics, channel_noise)
from test_core import synthetic


def props(sign, snr, channels):
    return dict(id='x', sign=sign, peak_snr=snr, channels=channels, pixels=10)


def test_narrow_component_like_off_source_artefacts_is_unreliable():
    # Five narrow (2-channel) emission artefacts far from the source, one narrow
    # and one broad emission component inside it.
    p = [props('Emission', 7., 2) for _ in range(5)] + [props('Emission', 6., 2), props('Emission', 5.5, 6)]
    radius = [2.]*5 + [.05, .05]
    areas = dict(inside=100, near=800, far=1000)
    narrow, broad = detection_reliability(p, radius, .1, areas)[-2:]
    assert narrow['reliability'] < .5 and narrow['control_far'] == 5
    assert broad['reliability'] > narrow['reliability'] and broad['control_far'] == 0
    assert detection_reliability(p, radius, .1, areas)[0]['region'] == 'Control'


def test_near_zone_artefacts_are_not_diluted_by_the_far_field():
    p = [props('Absorption', 9., 3)]*3 + [props('Absorption', 8., 3)]
    radius = [.2, .2, .2, .05]
    rel = detection_reliability(p, radius, .1, dict(inside=100, near=300, far=100000))[-1]
    # Near zone: 3 artefacts in 300 px -> 1 expected in the 100 px source region.
    assert rel['expected_false'] == pytest.approx(1.) and rel['reliability'] == 0.


def test_absorption_strength_is_positive_and_noise_is_robust():
    cube = synthetic()
    noise = channel_noise(cube)
    assert np.allclose(noise, .01, rtol=.1)
    comp = dict(id='A', sign='Absorption', lo=29, hi=33, x=40., y=40., footprint=[40*81+40, 40*81+41])
    p = component_properties(cube, [comp], noise)[0]
    assert p['peak_snr'] > 50 and p['channels'] == 5 and p['pixels'] == 2


def test_neighbours_within_one_beam():
    assert neighbour_counts([(0, 0), (.01, 0), (1, 1)], .02) == [1, 1, 0]


def test_identification_probabilities_and_chance_rate():
    (p,), none = identification_probabilities([0.], 30., window_kms=200.)
    assert p == pytest.approx(.842, abs=.002) and p + none == pytest.approx(1.)
    split, _ = identification_probabilities([0., 0.], 30.)
    assert split[0] == pytest.approx(split[1])
    (far,), none = identification_probabilities([120.], 30.)
    assert far < .01 and none > .99
    (shifted,), _ = identification_probabilities([40.], 30., mean_kms=40.)
    assert shifted == pytest.approx(p)
    # IC 860 SPW0 Extragalactic density: about 11 lines per GHz at 169 GHz.
    assert chance_match_probability(23, 2.094, 169., 100., k=2.) == pytest.approx(.92, abs=.01)
    assert chance_match_probability(23, 2.094, 169., 30.) == pytest.approx(.52, abs=.01)


def test_anchor_statistics_need_two_lines():
    assert anchor_statistics([10.]) is None
    stats = anchor_statistics([10., 12., 8., 30.])
    assert stats['n'] == 4 and stats['mean_kms'] == pytest.approx(11.)
