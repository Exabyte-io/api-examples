import pytest
from mat3ra.notebooks_utils.workflow import set_content

PW_INPUT = (
    "&system\n    VDW_CORR = 'dft-d'\n    ibrav = {{ input.IBRAV }}\n/\n"
    "&ELECTRONS\n    diago_full_acc = .false.\n    mixing_beta = 0.3\n/\n"
)


@pytest.mark.parametrize(
    "content,steps,present,absent,error",
    [
        (
            PW_INPUT,
            [
                ("system", {"vdw_corr": "d3_grimme"}),
                ("electrons", {"mixing_beta": 0.5, "diago_full_acc": True}),
            ],
            [
                "vdw_corr = 'd3_grimme'",
                "ibrav = {{ input.IBRAV }}",
                "mixing_beta = 0.5",
                "diago_full_acc = .true.",
            ],
            ["VDW_CORR = 'dft-d'", "mixing_beta = 0.3"],
            None,
        ),
        ("&SYSTEM\n/\n", [("IONS", {"ion_dynamics": "bfgs"})], [], [], "Namelist '&IONS' not found"),
    ],
)
def test_set_content(content, steps, present, absent, error):
    if error:
        with pytest.raises(ValueError, match=error):
            set_content(content, *steps[0])
        return
    result = content
    for section, parameters in steps:
        result = set_content(result, section, parameters)
    for text in present:
        assert text in result
    for text in absent:
        assert text not in result
