from pathlib import Path

import abangle.analyse as analyse


def test_analyse_uses_repo_data_directory():
    expected = Path(__file__).resolve().parents[1] / "data"

    assert Path(analyse.data_path) == expected


def test_analyse_uses_repo_config_fallback():
    expected = Path(__file__).resolve().parents[1] / "config"

    assert Path(analyse.config_path) == expected
