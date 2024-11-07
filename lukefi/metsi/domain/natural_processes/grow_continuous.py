from lukefi.metsi.forestry.naturalprocess.grow_continuous import continuous_growth_r
from lukefi.metsi.data.enums.internal import TreeSpecies
from lukefi.metsi.data.model import ForestStand, ReferenceTree


def grow_continuous(input: tuple[ForestStand, None], **operation_parameters) -> tuple[ForestStand, None]:
    step = operation_parameters.get('step', 5)
    stand, _ = input
    if len(stand.reference_trees) == 0:
        return input
    stand2 = continuous_growth_r(stand, step)
    return stand2, None
