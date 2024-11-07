from collections import defaultdict
import math
from statistics import median
from lukefi.metsi.data.model import ForestStand, ReferenceTree, TreeSpecies, SiteType, SoilPeatlandCategory

from lukefi.metsi.forestry.preprocessing.naslund import naslund_height


def yearly_diameter_growth_by_species(
    spe: TreeSpecies,
    d: float,
    h: float,
    biological_age_aggregate: float,
    d13_aggregate: float,
    height_aggregate: float,
    dominant_height: float,
    basal_area_total: float
) -> float:
    """ Model source: Acta Forestalia Fennica 163 """
    if spe == TreeSpecies.PINE:
        growth_percent = math.exp(5.4625
            - 0.6675 * math.log(biological_age_aggregate)
            - 0.4758 * math.log(basal_area_total)
            + 0.1173 * math.log(d13_aggregate)
            - 0.9442 * math.log(dominant_height)
            - 0.3631 * math.log(d)
            + 0.7762 * math.log(h))
    else:
        growth_percent = math.exp(6.9342
            - 0.8808 * math.log(biological_age_aggregate)
            - 0.4982 * math.log(basal_area_total)
            + 0.4159 * math.log(d13_aggregate)
            - 0.3865 * math.log(height_aggregate)
            - 0.6267 * math.log(d)
            + 0.1287 * math.log(h))
    return growth_percent


def yearly_height_growth_by_species(
    spe: TreeSpecies,
    d: float,
    h: float,
    biological_age_aggregate: float,
    d13_aggregate: float,
    height_aggregate: float,
    basal_area_total: float
) -> float:
    """ Model source: Acta Forestalia Fennica 163 """
    if spe == TreeSpecies.PINE:
        growth_percent = math.exp(5.4636
            - 0.9002 * math.log(biological_age_aggregate)
            + 0.5475 * math.log(d13_aggregate)
            - 1.1339 * math.log(h))
    else:
        growth_percent = (12.7402
            - 1.1786 * math.log(biological_age_aggregate)
            - 0.0937 * math.log(basal_area_total)
            - 0.1434 * math.log(d13_aggregate)
            - 0.8070 * math.log(height_aggregate)
            + 0.7563 * math.log(d)
            - 2.0522 * math.log(h))
    return growth_percent


def Pukkala_diameter_growth_by_species(
    spe: TreeSpecies,
    initialDiameter: float,
    G_plot: float,
    BAL_Total: float,
    BAL_Spruce: float,
    BAL_S_B: float,
    TS: int,
    sitetype: SiteType,
    soilpeat: SoilPeatlandCategory
) -> float:
# Parameters (Pukkala et al. 2021):                                           
    D_param = { 'Intercept':[-7.1552, -12.7527, -8.6306], 'sqrt_d':[0.4415, 0.1693, 0.5097], 'd':[-0.0685, -0.0301, -0.0829], \
            'ln_G_1':[-0.2027, -0.1875, -0.3864], 'BAL_Total':[-0.1236, -0.0563, 0], 'BAL_Spruce':[0, -0.0870, 0], \
            'BAL_Spruce_Broadleaf':[0, 0, -0.0545], 'ln_TS':[1.1198, 1.9747, 1.3163], 'Peat':[-0.2425, 0, 0], 'd_Pendula_or_Aspen':[0, 0, 0.0253], \
            'Fertility':{'Herb-rich':[0.1438,0.2688,0.2566], 'Mesic':[0,0,0], 'Sub-xeric':[-0.1754,-0.2145,-0.2256], 'Xeric':[-0.5163,-0.6179,-0.3237]} }

# re-codings
    if spe == TreeSpecies.PINE:
        spi = 0
    elif spe == TreeSpecies.SPRUCE:
        spi = 1
    else:
        spi = 2

    if sitetype == (SiteType.VERY_RICH_SITE or SiteType.RICH_SITE):
        site = 'Herb-rich'
    elif sitetype == SiteType.DAMP_SITE:
        site = 'Mesic'
    elif sitetype == SiteType.SUB_DRY_SITE:
        site = 'Sub-xeric'
    else:
        site = 'Xeric'

    if soilpeat != SoilPeatlandCategory.MINERAL_SOIL:
        peat = 1
    else:
        peat = 0

    ln_D_increment = (D_param['Intercept'][spi]
             + D_param['sqrt_d'][spi] * math.sqrt(initialDiameter)
             + D_param['d'][spi] * initialDiameter
             + D_param['ln_G_1'][spi] * math.log(G_plot+1)
             + D_param['BAL_Total'][spi] * (BAL_Total/math.sqrt(initialDiameter+1))
             + D_param['BAL_Spruce'][spi] * (BAL_Spruce/math.sqrt(initialDiameter+1)) 
             + D_param['BAL_Spruce_Broadleaf'][spi] * (BAL_S_B/math.sqrt(initialDiameter+1)) 
             + D_param['ln_TS'][spi] * math.log(TS) 
             + D_param['Fertility'][site][spi] 
             + D_param['Peat'][spi] * peat
             + D_param['d_Pendula_or_Aspen'][spi] * 0 * initialDiameter) # deciduous always assumed betula pubesc. (no information in forest data standard & small effect) 
    return math.exp(ln_D_increment)

def grow_diameter_and_height(
    stand: ForestStand,
    step: int = 5
) -> tuple[list[float], list[float]]:
    """ Diameter and height growth for trees with height > 1.3 meters. Based on Acta Forestalia Fennica 163. """
    if len(stand.reference_trees) == 0:
        return [], []
    trees = stand.reference_trees
#    group = defaultdict(list)
#    for i,t in enumerate(trees):
#        group[t.species].append(i)
    ds = [t.breast_height_diameter or 0 for t in trees]
    dspred = ds
    hs = [t.height for t in trees]
    hspred = hs
    sps = [t.species for t in trees]
    gs = [t.stems_per_ha*math.pi*(0.01*0.5*d)**2 for t,d in zip(trees,ds)]
    G = sum(gs)
    dspruces = [ds[jj] for jj,xx in enumerate(sps) if xx==TreeSpecies.SPRUCE]
    dsprucebroads = [ds[jj] for jj,xx in enumerate(sps) if xx != TreeSpecies.PINE]
    gspruces = [gs[jj] for jj,xx in enumerate(sps) if xx==TreeSpecies.SPRUCE]
    gsprucebroads = [gs[jj] for jj,xx in enumerate(sps) if xx != TreeSpecies.PINE]
    TS = 1322 # stand.degree_days attribute in metsi, but not included in forest data standard -> never inputted!
    sitetype = stand.site_type_category
    landclass = stand.soil_peatland_category

    for i,t in enumerate(trees):
        dcur = t.breast_height_diameter
        if dcur is None:
            dcur = 0
        bal = sum([gs[jj] for jj,xx in enumerate(ds) if xx>dcur]) # bal
        balspruces = sum([gspruces[jj] for jj,xx in enumerate(dspruces) if xx>dcur]) # bal | sp spruce
        balsprucebroads = sum([gsprucebroads[jj] for jj,xx in enumerate(dsprucebroads) if xx>dcur]) # bal | sp != pine

        if hs[i] >= 1.3:
            pd = Pukkala_diameter_growth_by_species(t.species, dcur, G, bal, balspruces, balsprucebroads, TS, sitetype, landclass)
            #pd = yearly_diameter_growth_by_species(spe, ds[i], hs[i], ag, dg, hg, hdom, G)/100
            #ph = yearly_height_growth_by_species(spe, ds[i], hs[i], ag, dg, hg, G)/100
            height1 = naslund_height(dcur, t.species)
            if height1 is None:
                height1 = 0
            height2 = naslund_height(dcur+pd, t.species)
            dspred[i] += pd
            #ds[i] *= 1+pd
            hspred[i] += height2-height1
            #hs[i] *= 1+ph
        else:
            if not dspred[i]:
                dspred[i] = 1.0
            if not hspred[i]:
                hspred[i] = 0.0
            coef = (5-(1.3-hspred[i])/0.3)/5 
            pd = coef * Pukkala_diameter_growth_by_species(t.species, dspred[i], G, bal, balspruces, balsprucebroads, TS, sitetype, landclass)
            height1 = naslund_height(dspred[i], t.species)
            height2 = naslund_height(dspred[i]+pd, t.species)
            dspred[i] += pd
            hspred[i] = min(hspred[i], 1.3)
            hspred[i] += height2-height1

    return dspred, hspred
