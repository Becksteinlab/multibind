import numpy as np


def protonation_free_energy(pKa: float, pH: float = 0) -> float:
    """Calculate the free energy to protonate a site.

    Parameters
    ----------
    pKa : float
        pKa of the interaction site
    pH : float, optional
        System pH. Default (0) gives the standard state protonation free energy.

    Returns
    -------
    float
        The protonation free energy of an interaction site.
    """
    return np.log(10) * (pH - pKa)


def protonation_free_energy_standard_error(pKa_error: float) -> float:
    """Standard error in the protonation free energy given the standard error of the pKa.

    Parameters
    ----------
    pKa_error : float
        The standard error of the pKa.

    Returns
    -------
    float
        The standard error in the protonation free energy.
    """
    return np.log(10) * pKa_error


def binding_free_energy_general(std_dG: float, concentration: float = 1, std_state: float = 1) -> float:
    """Binding free energy as a function of ligand concentration.

    Parameters
    ----------
    std_dG : float
        The standard state binding free energy of a ligand to an interaction site.
    concentration : float, optional
        The bulk concentration of the ligand in molar.
        The default of 1 M returns the standard state binding free energy.
    std_state: float, optional
        Standard state concentration. This is canonically 1 M.

    Returns
    -------
    float
        Binding free energy of a general ligand.
    """
    return std_dG - np.log(concentration / std_state)
