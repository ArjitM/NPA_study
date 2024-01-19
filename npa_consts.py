"""
Collection of various constants outlining names of NPA database columns by group and legend of one-hot coded checkbox
values.
"""
__author__ = "Arjit M; amisra2@illinois.edu"
__version__ = "Jan 14 2024"

"""
List of column names indicating whether specific follow-up symptoms have resolved.
Not all names have _res suffix.
"""
FUP_RES = [
    'fup_new_memory_res',
    'fup_headache_res',
    'fup_hydro_res',
    'fup_nausea_res',
    'fup_seizure_res',
    'fup_vai_res',
    'fup_vail_res',
    'fup_vair_res',
    'fup_vfl_res',
    'fup_vfll_res',
    'fup_vflr_res',
    'fup_vflb_res',
    'fup_dv_res',
    'fup_fn_res',
    'fup_fnl_res',
    'fup_fnr_res',
    'fup_fw_res',
    'fup_fwl_res',
    'fup_fwr_res',
    'fup_hearimp_res',
    'fup_dysphagia_res',
    'fup_vcd_res',
    'fup_tongue_res',
    'fup_wue_res',
    'fup_wuel_res',
    'fup_wuer_res',
    'fup_wle_res',
    'fup_wlel_res',
    'fup_wler_res',
    'fup_wtm_res',
    'fup_fnd_res',
    'fup_fsd_res',
    'fup_gait_res',
    'fup_other_res',
    'fup_dysphasia', #dysarthia
    'fup_galactorrhea',
    'fup_amenorrhea',
    'fup_vertigo',
    'fup_exper_aph',
    'fup_rec_aph',
    'fup_ams',
]


"""
List of column names indicating locations of specific follow-up symptoms.
All names have _res suffix.
"""
FUP_LOC = [
    'fup_vai_loc',
    'fup_vfl_loc',
    'fup_fw_loc',
    'fup_wue_loc',
    'fup_wle_loc',
    'fup_fn_loc',
]

"""
List of column names concerning tumor physical characteristics, size, etc.
"""
TUMOR_VARS = [
    'planes___1',
    'planes___2',
    'planes___3',
    'mid_shift',
    'tsize_diam1',
    'tsize_axial_unit_2',
    'tsize_orient1',
    'tsize_diam2',
    'tsize_axial_unit_4',
    'tsize_orient2',
    'tsize_diam3',
    'tsize_axial_unit_3',
    'tsize_orient3',
    'tsize_axial', # is max size
    'tsize_axial_unit',
    'tsize_coronal',
    'tsize_coronal_unit',
    'tsize_oblique',
    'tsize_oblique_unit',
    'tsize_sagittal',
    'tsize_sagittal_unit',
    'tumor_control',
    'tumor_laterality___1',
    'tumor_laterality___2',
    'tumor_laterality___3',
    'tumor_laterality___4',
    'tumor_laterality___5',
    'tumor_loc___1',
    'tumor_loc___10',
    'tumor_loc___11',
    'tumor_loc___12',
    'tumor_loc___13',
    'tumor_loc___14',
    'tumor_loc___15',
    'tumor_loc___16',
    'tumor_loc___17',
    'tumor_loc___18',
    'tumor_loc___19',
    'tumor_loc___2',
    'tumor_loc___3',
    'tumor_loc___4',
    'tumor_loc___5',
    'tumor_loc___6',
    'tumor_loc___7',
    'tumor_loc___8',
    'tumor_loc___9',
    'two_staged',
]


# Aggregate variables

"""
Components of the Promis-29 survey, grouped by aggregate output variable.
Each output variable is a aggregate of survey questions.
"""
P29_COMPONENTS = {}
P29_COMPONENTS["OUTPUTS"] = [
    'p29_pf_raw',
    'p29_anxiety_raw',
    'p29_depression_raw',
    'p29_fatigue_raw',
    'p29_sd_raw',
    'p29_social_raw',
    'p29_pain_raw',
]
P29_COMPONENTS["OUTPUTS_POS"] = [
    'p29_anxiety_raw',
    'p29_depression_raw',
    'p29_fatigue_raw',
    'p29_sd_raw',
    'p29_pain_raw',
]
P29_COMPONENTS["OUTPUTS_NEG"] = [
    'p29_pf_raw',
    'p29_social_raw',
]
P29_COMPONENTS["p29_pf_raw"] = [
    "p29_pfa11",
    "p29_pfa21",
    "p29_pfa23",
    "p29_pfa53",
]
P29_COMPONENTS["p29_anxiety_raw"] = [
    'p29_edanx01',
    'p29_edanx40',
    'p29_edanx41',
    'p29_edanx53',
]
P29_COMPONENTS["p29_depression_raw"] = [
    "p29_eddep04",
    "p29_eddep06",
    "p29_eddep29",
    "p29_eddep41",
]
P29_COMPONENTS["p29_fatigue_raw"] = [
    "p29_hi7",
    "p29_an3",
    "p29_fatexp41",
    "p29_fatexp40",
]
P29_COMPONENTS["p29_sd_raw"] = [
    "p29_sleep109",
    "p29_sleep116",
    "p29_sleep20",
    "p29_sleep44",
]
P29_COMPONENTS["p29_social_raw"] = [
    "p29_srpper11_caps",
    "p29_srpper18_caps",
    "p29_srpper23_caps",
    "p29_srpper46_caps",
]
P29_COMPONENTS["p29_pain_raw"] = [
    "p29_painin9",
    "p29_painin22",
    "p29_painin31",
    "p29_painin34",
]


# Legend for various labels 'KEYWORD___/d+'
# \s*\|\s*

APPROACH_LABELS = {
    1: "Frontal craniotomy",
    2: "Bifrontal craniotomy",
    3: "Supraorbital craniotomy",
    4: "Orbitozygomatic approach",
    5: "Orbitotomy",
    6: "Zygomatic osteotomy",
    7: "Pterional (frontotemporal) craniotomy",
    8: "Anterior clinoidectomy",
    9: "Middle fossa approach",
    10: "Retrosigmoid approach",
    11: "Endoscopic endonasal approach",
    12: "Expanded endoscopic endonasal approach",
    13: "Far lateral / Extreme lateral approach",
    14: "Midline suboccipital approach",
    15: "Para-sagittal convexity approach",
    16: "Anterior petrosectomy (Kawase Approach)",
    17: "Posterior petrosectomy / Posterior Petrosal approach",
    18: "Trans-labyrinthine approach",
    19: "Transmastoid approach",
    20: "Transotic approach",
    21: "Transcochlear approach",
    22: "Supracerebellar infratentorial approach",
    23: "Transcallosal transchoroidal approach",
    24: "Telovelar approach",
    25: "Other",
    26: "Parietal Craniotomy",
    27: "Occipital craniotomy",
    28: "Temporal craniotomy",
}
TUMOR_LAT_LABELS = {
    1: "Right",
    2: "Left",
    3: "Bilateral - one lesion crossing midline",
    4: "Bilateral - separate lesions",
    5: "Midline",
}
TUMOR_LOC_LABELS = {
    1: "Brainstem",
    2: "Cerebellum",
    3: "Corpus Callosum",
    4: "Frontal lobe",
    5: "Hypothalamus",
    6: "Limbic lobe",
    7: "Occipital lobe",
    8: "Parietal lobe",
    9: "Sellar region",
    10: "Temporal lobe",
    11: "Tentorial region",
    12: "Thalamus",
    13: "Lateral Ventricle",
    14: "Third Ventricle",
    15: "Fourth Ventricle",
    16: "Insula",
    17: "Cerebellopontine Angle",
    18: "Cranial Nerves (one or more)",
    19: "Pineal Region",
}
PREOP_TX_LABELS = {
    1: "Chemotherapy",
    2: "Embolization",
    3: "Hormonal therapy",
    4: "Immunotherapy",
    5: "Molecular targeted agent",
    6: "Radiation therapy (WBRT)",
    7: "Radiation therapy (IMRT)",
    8: "Radiation therapy (Other type)",
    9: "Stereotactic radiosurgery (SRS)",
    10: "Observation",
    11: "Steroids",
    12: "None",
}
INPATIENT_LABELS = {
    1: "Cerebral spinal fluid (CSF) leak",
    2: "Deep vein Thrombosis",
    3: "Myocardial infarction (MI)",
    4: "New neurological deficit",
    5: "Post-operative seizure < 14 days",
    6: "Post-operative ischemic stroke",
    7: "Post-operative hemorrhage",
    8: "Pulmonary Embolism",
    9: "Surgical site infection (SSI)",
    10: "Other",
    98: "None",
    99: "Information not available",
}
DISHARGE_DISP_LABELS = {
    1: "Home with no services",
    2: "Home with home healthcare services",
    3: "Inpatient rehab",
    4: "Post-acute or non-acute care (skilled nursing care)",
    5: "Transferred to another acute care facility (hospital)",
    6: "Against medical advice (AMA)",
    7: "Death",
    8: "Hospice/Palliative Care",
}
FUP_SYMPTOM_LABELS = {
    36: "Acromegaly",
    37: "Altered mental status",
    32: "Amenorrhea",
    10: "Double vision",
    30: "Dysarthria/Slurred speech",
    16: "Dysphagia / Trouble Swallowing",
    34: "Expressive aphasia",
    40: "Facial numbness",
    41: "Facial weakness",
    25: "Focal sensory deficit",
    26: "Gait impairment",
    31: "Galactorrhea",
    1: "Headache",
    15: "Hearing impairment",
    3: "Nausea/Dizziness",
    0: "New-onset memory loss",
    28: "No symptoms found on screening",
    27: "Other neurological deficit",
    35: "Receptive aphasia",
    4: "Seizure",
    2: "Symptoms attributed to Hydrocephalus",
    18: "Tongue weakness",
    24: "Tremor",
    99: "Unknown",
    33: "Vertigo",
    38: "Visual acuity impairment",
    39: "Visual field loss",
    17: "Vocal cord dysfunction",
    43: "Weakness of Lower extremity",
    23: "Weakness of trunk muscles",
    42: "Weakness of Upper extremity",
}
REASON_READMIT_LABELS = {
    1: "Cerebral spinal fluid (CSF) leak (including subcutaneous CSF collections",
    2: "Deep vein thrombosis",
    3: "Hydrocephalus",
    4: "Hemorrhage/Hematoma",
    5: "Medical reason",
    6: "Myocardial infarction (MI)",
    7: "New neurological deterioration",
    8: "Pulmonary embolism",
    9: "Seizure < 14 days postop",
    10: "Stroke/CVA",
    11: "Surgical site infection",
    12: "Tumor progression",
    13: "Other",
    99: "Information not available",
}

PREFIX_TO_LABELS = {
    "approach": APPROACH_LABELS,
    "tumor_laterality": TUMOR_LAT_LABELS,
    "tumor_loc": TUMOR_LOC_LABELS,
    "preoptx_type":PREOP_TX_LABELS,
    "inpatient_complications": INPATIENT_LABELS,
    "discharge_disp": DISHARGE_DISP_LABELS,
    "fup_symptoms": FUP_SYMPTOM_LABELS,
    "reason_readmit": REASON_READMIT_LABELS,
}

