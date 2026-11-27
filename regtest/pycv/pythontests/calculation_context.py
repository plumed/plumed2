import plumedCommunications as PLMD

plumedInit = {"Value": PLMD.defaults.COMPONENT_NODEV, "ATOMS": "1"}
functionInit = {"Value": PLMD.defaults.COMPONENT_NODEV}
counts = {}


def plumedCalculate(action):
    calculation_only = action.isCalculationOnly()
    if not calculation_only:
        counts[action.label] = counts.get(action.label, 0) + 1
    return (
        100 * counts.get(action.label, 0)
        + 10 * int(calculation_only)
        + int(action.isExchangeStep())
    )
