# HI (line) survey object

from __future__ import print_function

__author__ = "E.A.K. Adams"

"""
Survey object that collects survey information
Methods for changing spectral resolution and getting
updated survey parameters
"""

import astropy.units as u
from astropy.table import QTable
from astropy.io import ascii
import modules.functions as func
import numpy as np
import string


class Survey(object):
    """
    A class holding survey information
    
    Attributes
    ----------
    beam : Quantity
        Angular resolution
    chan_size : Quantity
        spectral resolution
    rms : Quantity
        noise for given spectral resolution
    name : string
        name of survey
    facility : string
        name of facililty
    redshift : string
        redshift range
    nhi : float
        Column density sensitivity

    Methods
    -------
    change_spec_res
    """

    def __init__(self, name, beam, chan_size, rms,
                 facility=None, sky_coverage=np.nan, target=None,
                 redshifts=None, references = None):
        """
        Construct initial attributes for Observations object

        Parameters
        ----------
        
        """
        self.beam = beam
        self.chan_size = chan_size
        self.rms = rms
        self.name = name
        self.facility = facility
        self.sky_coverage = sky_coverage
        self.target = target
        self.redshifts = redshifts
        self.references = references

        # calculate NHI from provided inputs
        self.nhi = func.get_nhi(self.rms, self.chan_size, self.beam, self.beam)

    def change_spec_res(self, new_res):
        """
        Change spectral resolution to new specifed resolution, updating all quantities
        """
        # update rms
        self.rms = func.get_sens_for_res(self.rms, self.chan_size, new_res)

        # update chan_size, only after have updated rms
        self.chan_size = new_res

        # update nhi
        self.nhi = func.get_nhi(self.rms, self.chan_size, self.beam, self.beam)


class Survey_Collection(object):
    """
    A class that collects multiple survey objects

    Attributes
    ----------
    survey_list : list
        List-(like) of Survey objects
    survey_table : Table
        astropoy table collecting all survey information
    common_spec_res : Quantity
        common spectral resolution to take all surveys to

    Methods
    -------
    set_common_spec_res
        Sets the common spec res across all survey objects
    write_latex_table
        Write a latex table of all the surveys at common spectral resolution
    """

    def __init__(self, survey_list, common_spec_res=25 * u.kHz):
        """
        Take all surveys to common_spec_res and
        Construct astropy table containing all survey info
        """
        # set common_spec_res as attribute
        self.common_spec_res = common_spec_res
        # set up row data
        survey_rows = []
        self.ref_dict = {}
        ref_letter_list = list(string.ascii_lowercase)
        try:
            for survey in survey_list:
                # change to common res
                survey.change_spec_res(common_spec_res)
                # handle refs
                if isinstance(survey.references, list):
                    for i, ref in enumerate(survey.references):
                        letter = ref_letter_list.pop(0)
                        self.ref_dict[letter] = ref
                        survey.references[i] = letter
                elif isinstance(survey.references, str):
                    letter = ref_letter_list.pop(0)
                    self.ref_dict[letter] = survey.references
                    survey.references = letter
                # append to survey list
                survey_rows.append((survey.name, survey.facility, survey.sky_coverage,
                                    survey.target, survey.beam,
                                    survey.rms, survey.nhi, survey.redshifts,
                                    survey.references))



        except:
            print("Could not parse list of survey objects for spec resolution change")

        try:
            self.survey_table = QTable(rows=survey_rows,
                                       names=['survey', 'facility', 'coverage', 'targets',
                                              'beam', 'rms', 'nhi', 'redshifts', 'references'])
        except:
            print("Could not make survey table table")

    def write_latex_table(self, latex_table='surveys.tex'):
        """
        Write out latex table to latex_table path

        Have this as a separate method so can customize and add inputs as needed
        """
        # and write the table out:
        col_names = ['Survey', 'Facility', 'Sky Coverage', 'Targets', 'Beam size', 'rms', '$N_{HI}$', 'Redshift', 'Refs']

        colalign = 'lllp{2.5cm}llll'
        self.survey_table['nhi_1e18'] = self.survey_table['nhi'] / 1e18

        tablefoot_text = r"\tablefoot{"
        for key, value in zip(self.ref_dict.keys(), self.ref_dict.values()):
            tablefoot_text += f"{key} : {value}, "
        tablefoot_text += "}"

        if 'Hz' in self.common_spec_res.unit.to_string():
            self.common_spec_res = func.convert_chan_freq_vel(self.common_spec_res)

        ascii.write(self.survey_table['survey', 'facility', 'coverage', 'targets', 'beam', 'rms', 'nhi_1e18',
                                      'redshifts', 'references'],
                    latex_table,
                    format='latex',
                    overwrite=True,
                    names=col_names,
                    col_align=colalign,
                    latexdict={'header_start': "\hline \hline",
                               'header_end': "\hline",
                               'data_end': "\hline",
                               'caption': f"HI surveys at common spectral resolution {self.common_spec_res}",
                               'tablefoot' : tablefoot_text,
                               'preamble': ["\centering", "\small", "\label{tab:hi_surveys}"]},
                    formats={'rms': '5.2f', '$N_{HI}$': '5.1f', 'Beam size': '3.0f', 'Sky Coverage': '5.0f'}
                    )
