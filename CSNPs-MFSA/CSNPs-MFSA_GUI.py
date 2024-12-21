"""
Last updated: 2024/08/23
CSNP Structure Annotation Tool
"""

import os
import re
import time
import threading
import pandas as pd
import sqlite3
import numpy as np
import tkinter as tk
import customtkinter as ctk

from typing import Dict, List, Tuple
from pyteomics import mass, mgf
from pyteomics.mass import Composition
from matchms import Spectrum, calculate_scores
from matchms.similarity import CosineGreedy
from tkinter import ttk, filedialog, messagebox
from concurrent.futures import ThreadPoolExecutor
from PIL import ImageTk
from rdkit import Chem
from rdkit.Chem import Draw

class MyTabView(ctk.CTkTabview):
    def __init__(self, master, **kwargs):
        super().__init__(master, **kwargs)

        self.specific_ion_sets = []

        # create tabs
        self.add("1.Composition Calculate")
        self.add("2.Target Extract")
        self.add("3.Adduct Ion Calculate")
        self.add("4.Neutral Loss Extract")
        self.add("5.Feature Ion Extract")
        self.add("6.Cosine Score Calculate")
        self.add("7.Annotation")
        self.add("8.Visualization")

        # add widgets on tabs
        self.create_mf_calculator_tab()
        self.create_extractor_tab()
        self.create_adduct_ion_calculator_tab()
        self.create_neutral_loss_extractor_tab()
        self.create_feature_ion_extractor_tab()
        self.create_cosine_score_calculator_tab()
        self.create_annotation_tab()
        self.create_visualization_tab()


    '''mf_calculator_tab'''
    def create_mf_calculator_tab(self):
        mf_calculator_tab = self.tab("1.Composition Calculate")

        # Create a canvas
        canvas = tk.Canvas(mf_calculator_tab)
        canvas.pack(side="left", fill="both", expand=True)

        # Create a frame inside the canvas
        frame_inside_canvas = ctk.CTkFrame(canvas)
        canvas.create_window((0, 0), window=frame_inside_canvas, anchor='nw')

        # CTkFrame of Import CSV
        frame_ImportCSV = ctk.CTkFrame(frame_inside_canvas, fg_color="transparent")
        frame_ImportCSV.pack(fill='x', padx=10, pady=5)

        label_csv = ctk.CTkLabel(frame_ImportCSV, text="Import query MS Spectra:", font=('Arial', 14, 'bold'))
        label_csv.pack(side='left', padx=(0, 10))
        self.entry_MF_cal_csv_path = ctk.CTkEntry(frame_ImportCSV, placeholder_text="Input CSV file path", width=500)
        self.entry_MF_cal_csv_path.pack(side='left', padx=(0, 5))

        self.button_MF_cal_csv_browse = ctk.CTkButton(frame_ImportCSV, text="Browse file", command=self.browse_mf_calculator_csv, fg_color="#6598D3", hover_color="#5578B3")
        self.button_MF_cal_csv_browse.pack(padx=10, pady=5, anchor='e')

        # CTkFrame of Import DB
        frame_ImportDB = ctk.CTkFrame(frame_inside_canvas, fg_color="transparent")
        frame_ImportDB.pack(fill='x', padx=10, pady=5)

        label_db = ctk.CTkLabel(frame_ImportDB, text="Import formula database (POS/NEG):", font=('Arial', 14, 'bold'))
        label_db.pack(side='left', padx=(0, 10))
        self.entry_MF_cal_db_path = ctk.CTkEntry(frame_ImportDB, placeholder_text="Input database file path", width=424)
        self.entry_MF_cal_db_path.pack(side='left', padx=(0, 5))

        self.button_MF_cal_db_browse = ctk.CTkButton(frame_ImportDB, text="Browse file", command=self.browse_db, fg_color="#6598D3", hover_color="#5578B3")
        self.button_MF_cal_db_browse.pack(padx=10, pady=5, anchor='e')

        # CTkFrame of Mass Tolerance (ppm) and Charge
        frame_MF_cal_charge_ppm = ctk.CTkFrame(frame_inside_canvas, fg_color="transparent")
        frame_MF_cal_charge_ppm.pack(fill='x', padx=10, pady=5)
        # Charge
        label_MF_cal_charge = ctk.CTkLabel(frame_MF_cal_charge_ppm, text="・Charge:", font=('Arial', 14))
        label_MF_cal_charge.pack(side='left', padx=(0, 10))
        self.entry_MF_cal_charge = ctk.CTkEntry(frame_MF_cal_charge_ppm, placeholder_text="Input charge (e.g., +1,-1)", width=275)
        self.entry_MF_cal_charge.pack(side='left', padx=(0, 20))
        # Mass Tolerance (ppm)
        label_MF_cal_ppm = ctk.CTkLabel(frame_MF_cal_charge_ppm, text="・Mass Tolerance (ppm):", font=('Arial', 14))
        label_MF_cal_ppm.pack(side='left', padx=(0, 10))
        self.entry_MF_cal_ppm = ctk.CTkEntry(frame_MF_cal_charge_ppm, placeholder_text="Input ppm (e.g., 5)", width=275)
        self.entry_MF_cal_ppm.pack(side='left')

        # CTkFrame of Output
        frame_MF_cal_Output = ctk.CTkFrame(frame_inside_canvas, fg_color="transparent")
        frame_MF_cal_Output.pack(fill='x', padx=10, pady=5)
        label_MF_cal_output = ctk.CTkLabel(frame_MF_cal_Output, text="Export Result Files:", font=('Arial', 14, 'bold'))
        label_MF_cal_output.pack(side='left', padx=(0, 10))
        self.entry_MF_cal_output_path = ctk.CTkEntry(frame_MF_cal_Output, placeholder_text="Output csv file path, e.g. C:/Users/xxx/", width=330)
        self.entry_MF_cal_output_path.pack(side='left', padx=(0, 5))
        self.entry_MF_cal_output_filename = ctk.CTkEntry(frame_MF_cal_Output, placeholder_text="Output csv file name, e.g. xxx.csv", width=205)
        self.entry_MF_cal_output_filename.pack(side='left', padx=(0, 5))
        self.button_MF_cal_calculate = ctk.CTkButton(frame_MF_cal_Output, text="Calculate", command=self.calculate, fg_color="#6598D3", hover_color="#5578B3")
        self.button_MF_cal_calculate.pack(padx=10, pady=5, anchor='e')

    def browse_mf_calculator_csv(self):
        file_path = filedialog.askopenfilename(filetypes=[("CSV files", "*.csv")])
        if file_path:
            self.entry_MF_cal_csv_path.delete(0, tk.END)
            self.entry_MF_cal_csv_path.insert(0, file_path)

    def browse_db(self):
        file_path = filedialog.askopenfilename(filetypes=[("Database files", "*.db")])
        if file_path:
            self.entry_MF_cal_db_path.delete(0, tk.END)
            self.entry_MF_cal_db_path.insert(0, file_path)

    def calculate(self):
        csv_path = self.entry_MF_cal_csv_path.get()
        db_path = self.entry_MF_cal_db_path.get()
        output_path = self.entry_MF_cal_output_path.get()
        output_filename = self.entry_MF_cal_output_filename.get()

        if os.path.exists(csv_path) and os.path.exists(db_path) and output_path and output_filename:
            df, output_file = self.process_mf_calculator(csv_path, db_path, output_path, output_filename)
            df.to_csv(output_file, index=False)
            messagebox.showinfo(title="Success", message=f"Data processed and saved to {output_file}")
        else:
            messagebox.showerror(title="Error", message="Please make sure all file paths and names are correct and accessible.")

    def process_mf_calculator(self, csv_path, db_path, output_path, output_filename):
        charge = int(self.entry_MF_cal_charge.get())
        ppm = float(self.entry_MF_cal_ppm.get())
        formula_ranges = self.read_MF_cal_formulas_from_db(db_path, charge, ppm)
        df = pd.read_csv(csv_path)
        elemental_compositions = df['row m/z'].apply(lambda mz: self.find_MF_cal_matching_formulas(mz, formula_ranges))
        mz_col_index = df.columns.get_loc('row m/z') + 1
        df.insert(mz_col_index, 'elemental composition', elemental_compositions)
        output_csv_path = os.path.join(output_path, output_filename)
        return df, output_csv_path

    def mf_calculator_mz_range(self, formula: str, charge: int, ppm: float) -> tuple:
        mz = mass.calculate_mass(formula=formula, charge=charge)
        mz_min = mz - mz * ppm / 1e6
        mz_max = mz + mz * ppm / 1e6
        return abs(mz_min), abs(mz_max)

    def read_MF_cal_formulas_from_db(self, db_path: str, charge: int, ppm: float) -> Dict[str, Tuple[float, float]]:
        conn = sqlite3.connect(db_path)
        cursor = conn.cursor()
        cursor.execute("SELECT formula FROM compounds")
        formulas = cursor.fetchall()
        conn.close()
        return {formula[0]: self.mf_calculator_mz_range(formula[0], charge, ppm) for formula in formulas}

    def find_MF_cal_matching_formulas(self, mz_value: float, formula_ranges: dict) -> str:
        min_ppm = float('inf')
        best_match = ''
        for formula, (mz_min, mz_max) in formula_ranges.items():
            if mz_min <= mz_value <= mz_max:
                mz = (mz_max + mz_min) / 2
                ppm = abs((mz_value - mz) / mz * 1e6)
                if ppm < min_ppm:
                    min_ppm = ppm
                    best_match = formula
        return best_match


    '''extractor_tab'''
    def create_extractor_tab(self):
        extractor_tab = self.tab("2.Target Extract")

        # Create a canvas
        canvas = tk.Canvas(extractor_tab)
        canvas.pack(side="left", fill="both", expand=True)

        # Create a frame inside the canvas
        frame_inside_canvas = ctk.CTkFrame(canvas)
        canvas.create_window((0, 0), window=frame_inside_canvas, anchor='nw')

        # CTkFrame of Import MGF
        frame_ImportMGF = ctk.CTkFrame(frame_inside_canvas, fg_color="transparent")
        frame_ImportMGF.pack(fill='x', padx=10, pady=5)
        # Section for MGF file input
        label_extractor_mgf = ctk.CTkLabel(frame_ImportMGF, text="Import MGF file:", font=('Arial', 14, 'bold'))
        label_extractor_mgf.pack(side='left', padx=(0, 10))
        self.entry_extractor_mgf_path = ctk.CTkEntry(frame_ImportMGF, placeholder_text="Input file path", width=550)
        self.entry_extractor_mgf_path.pack(side='left', padx=(0, 5))
        self.button_extractor_browse = ctk.CTkButton(frame_ImportMGF, text="Browse file",
                                                    command=self.browse_extractor_mgf, fg_color="#6598D3", hover_color="#5578B3")
        self.button_extractor_browse.pack(padx=10, pady=5, anchor='e')

        # CTkFrame of Import CSV
        frame_ImportCSV = ctk.CTkFrame(frame_inside_canvas, fg_color="transparent")
        frame_ImportCSV.pack(fill='x', padx=10, pady=5)
        # Section for CSV file input
        label_extractor_csv = ctk.CTkLabel(frame_ImportCSV, text="Import CSV file:", font=('Arial', 14, 'bold'))
        label_extractor_csv.pack(side='left', padx=(0, 10))
        self.entry_extractor_csv_path = ctk.CTkEntry(frame_ImportCSV, placeholder_text="Input CSV file path", width=550)
        self.entry_extractor_csv_path.pack(side='left', padx=(0, 5))
        self.button_csv_browse = ctk.CTkButton(frame_ImportCSV, text="Browse file",
                                            command=self.browse_extractor_csv, fg_color="#6598D3", hover_color="#5578B3")
        self.button_csv_browse.pack(padx=10, pady=5, anchor='e')

        # TypeI_feature_formulas
        frame_TypeI_feature_formulas = ctk.CTkFrame(frame_inside_canvas, fg_color="transparent")
        frame_TypeI_feature_formulas.pack(fill='x', padx=10, pady=5)
        label_TypeI_feature = ctk.CTkLabel(frame_TypeI_feature_formulas, text="Type I Feature Formulas:", font=('Arial', 14, 'bold'))
        label_TypeI_feature.pack(side='left', padx=(0, 10))
        self.entry_TypeI_feature_formulas = ctk.CTkEntry(frame_TypeI_feature_formulas, placeholder_text="Input Type I feature ion formulas, e.g.,Daphnane, C20H28O4, C20H28O5", width=660)
        self.entry_TypeI_feature_formulas.pack(side='left', padx=(0, 5))

        # TypeI_feature_formulas; mass range of product ions
        frame_typeI_para1 = ctk.CTkFrame(frame_inside_canvas, fg_color="transparent")
        frame_typeI_para1.pack(fill='x', padx=10, pady=5)

        label_massrange_min  = ctk.CTkLabel(frame_typeI_para1, text="・Mass Range of Product Ions (m/z):", font=('Arial', 14))
        label_massrange_min.pack(side='left', padx=(0, 10))
        self.entry_mass_range_min_typeI = ctk.CTkEntry(frame_typeI_para1, placeholder_text="min, e.g.230", width=203)
        self.entry_mass_range_min_typeI.pack(side='left', padx=(0, 5))

        label_massrange_max  = ctk.CTkLabel(frame_typeI_para1, text="to", font=('Arial', 14))
        label_massrange_max.pack(side='left', padx=(0, 5))

        # Option Menu for max_mass
        self.max_mass_option = ctk.CTkOptionMenu(frame_typeI_para1, values=["Use Precursor Ion", "Specify Value"],
                                                fg_color="#F5F5F5",  # Background color
                                                button_color="#BEBEBE",  # Button color
                                                text_color="#333333",  # Text color
                                                width=160)
        self.max_mass_option.pack(side='left', padx=(0, 5))
        self.entry_mass_range_max_typeI = ctk.CTkEntry(frame_typeI_para1, placeholder_text="max, e.g.400", width=203)
        self.entry_mass_range_max_typeI.pack(side='left', padx=(0, 5))

        # TypeI_feature_formulas; mass range of product ions
        frame_typeI_para2 = ctk.CTkFrame(frame_inside_canvas, fg_color="transparent")
        frame_typeI_para2.pack(fill='x', padx=10, pady=5)

        label_charge_typeI  = ctk.CTkLabel(frame_typeI_para2, text="・Charge:", font=('Arial', 14))
        label_charge_typeI.pack(side='left', padx=(0, 5))
        self.entry_charge_typeI = ctk.CTkEntry(frame_typeI_para2, placeholder_text="eg.+1/-1", width=80)
        self.entry_charge_typeI.pack(side='left', padx=(0, 5))

        label_massTolerance_typeI  = ctk.CTkLabel(frame_typeI_para2, text="・Mass Tolerence (ppm):", font=('Arial', 14))
        label_massTolerance_typeI.pack(side='left', padx=(0, 5))
        self.entry_mass_tolerance_typeI = ctk.CTkEntry(frame_typeI_para2, placeholder_text="e.g., 5", width=80)
        self.entry_mass_tolerance_typeI.pack(side='left', padx=(0, 5))

        label_hitScore_typeI  = ctk.CTkLabel(frame_typeI_para2, text="・Hit Score:", font=('Arial', 14))
        label_hitScore_typeI.pack(side='left', padx=(0, 5))
        self.entry_hitScore_typeI = ctk.CTkEntry(frame_typeI_para2, placeholder_text="e.g., 5", width=80)
        self.entry_hitScore_typeI.pack(side='left', padx=(0, 5))

        label_name_typeI  = ctk.CTkLabel(frame_typeI_para2, text="・Type Name:", font=('Arial', 14))
        label_name_typeI.pack(side='left', padx=(0, 5))
        self.entry_name_typeI = ctk.CTkEntry(frame_typeI_para2, placeholder_text="e.g., Daphnane", width=180)
        self.entry_name_typeI.pack(side='left', padx=(0, 5))

        # TypeII_feature_formulas
        frame_TypeII_feature_formulas = ctk.CTkFrame(frame_inside_canvas, fg_color="transparent")
        frame_TypeII_feature_formulas.pack(fill='x', padx=10, pady=5)
        label_TypeII_feature = ctk.CTkLabel(frame_TypeII_feature_formulas, text="Type II Feature Formulas:", font=('Arial', 14, 'bold'))
        label_TypeII_feature.pack(side='left', padx=(0, 10))
        self.entry_TypeII_feature_formulas = ctk.CTkEntry(frame_TypeII_feature_formulas, placeholder_text="Input Type II feature ion formulas, e.g., Macrocyclic Daphnane, C29H42O8, C29H42O9", width=660)
        self.entry_TypeII_feature_formulas.pack(side='left', padx=(0, 5))

        # TypeII_feature_formulas; mass range of product ions
        frame_typeII_para1 = ctk.CTkFrame(frame_inside_canvas, fg_color="transparent")
        frame_typeII_para1.pack(fill='x', padx=10, pady=5)

        label_massrange_min  = ctk.CTkLabel(frame_typeII_para1, text="・Mass Range of Product Ions (m/z):", font=('Arial', 14))
        label_massrange_min.pack(side='left', padx=(0, 10))
        self.entry_mass_range_min_typeII = ctk.CTkEntry(frame_typeII_para1, placeholder_text="min, e.g., 130", width=203)
        self.entry_mass_range_min_typeII.pack(side='left', padx=(0, 5))

        label_massrange_max  = ctk.CTkLabel(frame_typeII_para1, text="to", font=('Arial', 14))
        label_massrange_max.pack(side='left', padx=(0, 5))

        # Option Menu for max_mass
        self.max_mass_option = ctk.CTkOptionMenu(frame_typeII_para1, values=["Use Precursor Ion", "Specify Value"],
                                                fg_color="#F5F5F5",  # Background color
                                                button_color="#BEBEBE",  # Button color
                                                text_color="#333333",  # Text color
                                                width=160)
        self.max_mass_option.pack(side='left', padx=(0, 5))

        self.entry_mass_range_max_typeII = ctk.CTkEntry(frame_typeII_para1, placeholder_text="max, e.g.400", width=203)
        self.entry_mass_range_max_typeII.pack(side='left', padx=(0, 5))

        # TypeII_feature_formulas; mass range of product ions
        frame_typeII_para2 = ctk.CTkFrame(frame_inside_canvas, fg_color="transparent")
        frame_typeII_para2.pack(fill='x', padx=10, pady=5)

        label_charge_typeII  = ctk.CTkLabel(frame_typeII_para2, text="・Charge:", font=('Arial', 14))
        label_charge_typeII.pack(side='left', padx=(0, 5))
        self.entry_charge_typeII = ctk.CTkEntry(frame_typeII_para2, placeholder_text="e.g., +1/-1", width=80)
        self.entry_charge_typeII.pack(side='left', padx=(0, 5))

        label_massTolerance_typeII  = ctk.CTkLabel(frame_typeII_para2, text="・Mass Tolerence (ppm):", font=('Arial', 14))
        label_massTolerance_typeII.pack(side='left', padx=(0, 5))
        self.entry_mass_tolerance_typeII = ctk.CTkEntry(frame_typeII_para2, placeholder_text="e.g., 5", width=80)
        self.entry_mass_tolerance_typeII.pack(side='left', padx=(0, 5))

        label_hitScore_typeII  = ctk.CTkLabel(frame_typeII_para2, text="・Hit Score:", font=('Arial', 14))
        label_hitScore_typeII.pack(side='left', padx=(0, 5))
        self.entry_hitScore_typeII = ctk.CTkEntry(frame_typeII_para2, placeholder_text="e.g.,4", width=80)
        self.entry_hitScore_typeII.pack(side='left', padx=(0, 5))

        label_name_typeII  = ctk.CTkLabel(frame_typeII_para2, text="・Type Name:", font=('Arial', 14))
        label_name_typeII.pack(side='left', padx=(0, 5))
        self.entry_name_typeII = ctk.CTkEntry(frame_typeII_para2, placeholder_text="e.g., Macrocyclic Daphnane", width=180)
        self.entry_name_typeII.pack(side='left', padx=(0, 5))

        # CTkFrame of Output1
        frame_Output1 = ctk.CTkFrame(frame_inside_canvas, fg_color="transparent")
        frame_Output1.pack(fill='x', padx=10, pady=5)
        label_output = ctk.CTkLabel(frame_Output1, text="Export Result Files:", font=('Arial', 14, 'bold'))
        label_output.pack(side='left', padx=(0, 10))
        self.entry_extractor_output_path = ctk.CTkEntry(frame_Output1,
                                                        placeholder_text="Output directory path, e.g. C:/Users/xxx/", width=680)
        self.entry_extractor_output_path.pack(side='left', padx=(0, 5))

        # CTkFrame of Output2
        frame_Output2 = ctk.CTkFrame(frame_inside_canvas, fg_color="transparent")
        frame_Output2.pack(fill='x', padx=10, pady=5)
        self.entry_extractor_output_mgf_filename = ctk.CTkEntry(frame_Output2,
                                                                placeholder_text="Output MGF file name, e.g. result.mgf", width=335)
        self.entry_extractor_output_mgf_filename.pack(side='left', padx=(0, 5))
        self.entry_output_csv_filename = ctk.CTkEntry(frame_Output2,
                                                    placeholder_text="Output CSV file name, e.g. result.csv", width=335)
        self.entry_output_csv_filename.pack(side='left', padx=(0, 5))
        self.button_extractor_calculate = ctk.CTkButton(frame_Output2, text="Extract and Classify",
                                                        command=self.extract_and_classify, fg_color="#6598D3", hover_color="#5578B3")
        self.button_extractor_calculate.pack(padx=10, pady=5, anchor='e')

        # Progress Bar
        # self.progress_bar = ctk.CTkProgressBar(frame_Output2)
        # self.progress_bar.pack(fill='x', padx=10)
        # self.progress_bar.set(0)

    def browse_extractor_mgf(self):
        file_path = filedialog.askopenfilename(filetypes=[("MGF files", "*.mgf")])
        if file_path:
            self.entry_extractor_mgf_path.delete(0, tk.END)
            self.entry_extractor_mgf_path.insert(0, file_path)

    def browse_extractor_csv(self):
        file_path = filedialog.askopenfilename(filetypes=[("CSV files", "*.csv")])
        if file_path:
            self.entry_extractor_csv_path.delete(0, tk.END)
            self.entry_extractor_csv_path.insert(0, file_path)

    def extract_and_classify(self):
        threading.Thread(target=self.run_extraction_and_classification).start()

    def run_extraction_and_classification(self):
        mgf_path = self.entry_extractor_mgf_path.get()
        csv_path = self.entry_extractor_csv_path.get()
        output_path = self.entry_extractor_output_path.get()
        output_mgf_filename = self.entry_extractor_output_mgf_filename.get()
        output_csv_filename = self.entry_output_csv_filename.get()
        if os.path.exists(mgf_path) and os.path.exists(
                csv_path) and output_path and output_mgf_filename and output_csv_filename:
            output_mgf_file = os.path.join(output_path, output_mgf_filename)
            output_csv_file = os.path.join(output_path, output_csv_filename)
            self.extract_and_classify_features(mgf_path, csv_path, output_mgf_file, output_csv_file)
            self.progress_bar.set(100)
            messagebox.showinfo(title="Success",
                                message=f"Data processed and saved to {output_mgf_file} and {output_csv_file}")
        else:
            messagebox.showerror(title="Error",
                                message="Please make sure all file paths and names are correct and accessible.")

        # Update progress bar in steps
        # for i in range(100):
            # time.sleep(0.1)  # Simulate processing time
            # self.progress_bar.set(i + 1)  # Update progress bar
            # self.update()  # Important: Update the GUI

    def extractor_mz_range(self, formula: str, *, charge: int, tolerance_ppm: float) -> tuple[float, float]:
        mz = mass.calculate_mass(formula=formula, charge=charge)
        mz_min = mz - mz * tolerance_ppm / 1e6
        mz_max = mz + mz * tolerance_ppm / 1e6
        return mz_min, mz_max

    def detect_TypeI_spectra(self, filename: str, feature_formulas: List[str]) -> Tuple[List[dict], List[int]]:
        try:
            min_mass = float(self.entry_mass_range_min_typeI.get())
            charge = int(self.entry_charge_typeI.get())
            tolerance_ppm = float(self.entry_mass_tolerance_typeI.get())
            max_mass_option = self.max_mass_option.get()
            if max_mass_option == "Specify Value":
                max_mass = float(self.entry_mass_range_max_typeI.get())
        except ValueError:
            messagebox.showerror("Input Error",
                                "Please enter valid numeric values for mass ranges, charge, and tolerance.")
            return [], []

        TypeI_spectra = []
        TypeI_HitScores_list = []
        with mgf.read(filename) as spectra:
            for spectrum in spectra:
                if max_mass_option == "Use Precursor Ion":
                    max_mass = spectrum['params']['pepmass'][0]
                product_ions = spectrum['m/z array']
                intensity = spectrum['intensity array']
                ind = np.nonzero((product_ions >= min_mass) & (product_ions <= max_mass))[0]
                if len(ind) == 0:
                    continue
                intensity_max = np.max(intensity[ind])
                intensity_min = np.min(intensity[ind])
                x = (intensity[ind] - intensity_min) / (intensity_max - intensity_min)
                Nor_intensity = np.zeros(intensity.shape[0])
                Nor_intensity[ind] = x
                TypeI_HitScores = sum(
                    1 for intensity_value, product_ion in zip(Nor_intensity, product_ions)
                    if intensity_value > 0.1 and any(
                        self.extractor_mz_range(formula, charge=charge, tolerance_ppm=tolerance_ppm)[0]
                        <= product_ion <=
                        self.extractor_mz_range(formula, charge=charge, tolerance_ppm=tolerance_ppm)[1]
                        for formula in feature_formulas)
                )
                if TypeI_HitScores >= int(self.entry_hitScore_typeI.get()):
                    TypeI_spectra.append(spectrum)
                    TypeI_HitScores_list.append(TypeI_HitScores)
        return TypeI_spectra, TypeI_HitScores_list

    def detect_TypeII_spectra(self, filename: str, feature_formulas: List[str]) -> Tuple[List[dict], List[int]]:
        try:
            min_mass = float(self.entry_mass_range_min_typeII.get())
            charge = int(self.entry_charge_typeII.get())
            tolerance_ppm = float(self.entry_mass_tolerance_typeII.get())
            max_mass_option = self.max_mass_option.get()
            if max_mass_option == "Specify Value":
                max_mass = float(self.entry_mass_range_max_typeII.get())
        except ValueError:
            messagebox.showerror("Input Error",
                                "Please enter valid numeric values for mass ranges, charge, and tolerance.")
            return [], []

        TypeII_spectra = []
        TypeII_HitScores_list = []
        with mgf.read(filename) as spectra:
            for spectrum in spectra:
                if max_mass_option == "Use Precursor Ion":
                    max_mass = spectrum['params']['pepmass'][0]
                product_ions = spectrum['m/z array']
                intensity = spectrum['intensity array']
                ind = np.nonzero((product_ions >= min_mass) & (product_ions <= max_mass))[0]
                if len(ind) == 0:
                    continue
                intensity_max = np.max(intensity[ind])
                intensity_min = np.min(intensity[ind])
                x = (intensity[ind] - intensity_min) / (intensity_max - intensity_min)
                Nor_intensity = np.zeros(intensity.shape[0])
                Nor_intensity[ind] = x
                TypeII_HitScores = sum(
                    1 for intensity_value, product_ion in zip(Nor_intensity, product_ions)
                    if intensity_value > 0.1 and any(
                        self.extractor_mz_range(formula, charge=charge, tolerance_ppm=tolerance_ppm)[0]
                        <= product_ion <=
                        self.extractor_mz_range(formula, charge=charge, tolerance_ppm=tolerance_ppm)[1]
                        for formula in feature_formulas)
                )
                if TypeII_HitScores >= int(self.entry_hitScore_typeI.get()):
                    TypeII_spectra.append(spectrum)
                    TypeII_HitScores_list.append(TypeII_HitScores)
        return TypeII_spectra, TypeII_HitScores_list

    def extract_and_classify_features(self, mgf_path, csv_path, output_mgf_file, output_csv_file):

        TypeI_featureMF = [x.strip() for x in self.entry_TypeI_feature_formulas.get().split(',')]
        TypeII_featureMF = [x.strip() for x in self.entry_TypeII_feature_formulas.get().split(',')]

        # TypeI
        TypeI_spectra, TypeI_HitScores = self.detect_TypeI_spectra(mgf_path, TypeI_featureMF)
        TypeI_IDs = [int(spectrum['params']['title']) for spectrum in TypeI_spectra]

        # TypeII
        TypeII_spectra, TypeII_HitScores = self.detect_TypeII_spectra(mgf_path, TypeII_featureMF)
        TypeII_IDs = [int(spectrum['params']['title']) for spectrum in TypeII_spectra]

        # Combine all spectra
        all_spectra = TypeI_spectra + TypeII_spectra

        # Remove duplicates
        unique_spectra = []
        seen_titles = set()
        for spectrum in all_spectra:
            title = spectrum['params']['title']
            if title not in seen_titles:
                seen_titles.add(title)
                unique_spectra.append(spectrum)

        # Sort by 'row ID'
        unique_spectra.sort(key=lambda x: int(x['params']['title']))

        # Write to a single mgf file
        mgf.write(unique_spectra, output=output_mgf_file)

        # Load a CSV file into a DataFrame with 'row ID' as the index.
        try:
            input_csv = pd.read_csv(csv_path, index_col='row ID')
        except pd.errors.ParserError:
            input_csv = pd.read_csv(csv_path, index_col='row ID', error_bad_lines=False)
        input_IDs = input_csv.index.tolist()
        drop_IDs = [x for x in input_IDs if x not in TypeI_IDs and x not in TypeII_IDs]
        focal_csv = input_csv.drop(drop_IDs, axis=0)
        focal_csv['TypeI Score'] = pd.Series({id_: score for id_, score in zip(TypeI_IDs, TypeI_HitScores)}, dtype=object).reindex(
            focal_csv.index)
        focal_csv['TypeII Score'] = pd.Series({id_: score for id_, score in zip(TypeII_IDs, TypeII_HitScores)}, dtype=object).reindex(
            focal_csv.index)
        name_typeI = self.entry_name_typeI.get()
        name_typeII = self.entry_name_typeII.get()
        focal_csv['type'] = np.where(focal_csv['TypeI Score'].isnull() | (focal_csv['TypeII Score'] > focal_csv['TypeI Score']), name_typeII, name_typeI)
        focal_csv = focal_csv.drop(['TypeI Score', 'TypeII Score'], axis=1)
        focal_csv = focal_csv.loc[:, ~focal_csv.columns.str.contains('^Unnamed')]
        focal_csv.to_csv(output_csv_file, encoding='UTF-8')


    '''adduct_ion_calculator_tab'''
    def create_adduct_ion_calculator_tab(self):
        adduct_ion_calculator_tab = self.tab("3.Adduct Ion Calculate")

        # Create a canvas
        canvas = tk.Canvas(adduct_ion_calculator_tab)
        canvas.pack(side="left", fill="both", expand=True)

        # Create a frame inside the canvas
        frame_inside_canvas = ctk.CTkFrame(canvas)
        canvas.create_window((0, 0), window=frame_inside_canvas, anchor='nw')

        # CTkFrame of Import CSV
        frame_ImportCSV = ctk.CTkFrame(frame_inside_canvas, fg_color="transparent")
        frame_ImportCSV.pack(fill='x', padx=10, pady=5)
        label_csv = ctk.CTkLabel(frame_ImportCSV, text="Import query MS Spectra:", font=('Arial', 14, 'bold'))
        label_csv.pack(side='left', padx=(0, 10))
        self.entry_adduct_csv_path = ctk.CTkEntry(frame_ImportCSV, placeholder_text="Input CSV file path", width=500)
        self.entry_adduct_csv_path.pack(side='left', padx=(0, 5))

        self.button_adduct_csv_browse = ctk.CTkButton(frame_ImportCSV, text="Browse file",
                                                    command=self.browse_adduct_csv, fg_color="#6598D3", hover_color="#5578B3")
        self.button_adduct_csv_browse.pack(padx=10, pady=5, anchor='e')

        # CTkFrame of Mass Tolerance (ppm) and Charge
        frame_adduct_charge_ppm = ctk.CTkFrame(frame_inside_canvas, fg_color="transparent")
        frame_adduct_charge_ppm.pack(fill='x', padx=10, pady=5)
        # Charge
        label_charge = ctk.CTkLabel(frame_adduct_charge_ppm, text="・Charge:", font=('Arial', 14))
        label_charge.pack(side='left', padx=(0, 10))
        self.entry_adduct_charge = ctk.CTkEntry(frame_adduct_charge_ppm, placeholder_text="Input charge (e.g., +1,-1)", width=278)
        self.entry_adduct_charge.pack(side='left', padx=(0, 20))
        # Mass Tolerance (ppm)
        label_ppm = ctk.CTkLabel(frame_adduct_charge_ppm, text="・Mass Tolerance (ppm):", font=('Arial', 14))
        label_ppm.pack(side='left', padx=(0, 10))
        self.entry_adduct_ppm = ctk.CTkEntry(frame_adduct_charge_ppm, placeholder_text="Input ppm (e.g., 5)", width=278)
        self.entry_adduct_ppm.pack(side='left')

        # Output File Path Section
        frame_adduct_Output = ctk.CTkFrame(frame_inside_canvas, fg_color="transparent")
        frame_adduct_Output.pack(fill='x', padx=10, pady=5)
        label_adduct_output = ctk.CTkLabel(frame_adduct_Output , text="Export Result Files:", font=('Arial', 14, 'bold'))
        label_adduct_output.pack(side='left', padx=(0, 10))
        self.entry_adduct_output_path = ctk.CTkEntry(frame_adduct_Output ,
                                                    placeholder_text="Output csv file path, e.g. C:/Users/xxx/", width=330)
        self.entry_adduct_output_path.pack(side='left', padx=(0, 5))
        self.entry_adduct_output_csv_filename = ctk.CTkEntry(frame_adduct_Output ,
                                                            placeholder_text="Output CSV file name, e.g. result.csv", width=205)
        self.entry_adduct_output_csv_filename.pack(side='left', padx=(0, 5))
        self.button_adduct_calculate = ctk.CTkButton(frame_adduct_Output , text="Calculate",
                                                    command=self.calculate_adduct, fg_color="#6598D3", hover_color="#5578B3")
        self.button_adduct_calculate.pack(padx=10, pady=5, anchor='e')

    def browse_adduct_csv(self):
        file_path = filedialog.askopenfilename(filetypes=[("CSV files", "*.csv")])
        if file_path:
            self.entry_adduct_csv_path.delete(0, tk.END)
            self.entry_adduct_csv_path.insert(0, file_path)

    def calculate_adduct(self):
        csv_path = self.entry_adduct_csv_path.get()
        output_path = self.entry_adduct_output_path.get()
        output_csv_filename = self.entry_adduct_output_csv_filename.get()
        charge = int(self.entry_adduct_charge.get())
        tolerance_ppm = float(self.entry_adduct_ppm.get())
        if os.path.exists(csv_path) and output_path:
            try:
                df = pd.read_csv(csv_path)
                df_processed = self.process_dataframe(df, charge, tolerance_ppm)
                output_file = os.path.join(output_path, output_csv_filename)
                df_processed.to_csv(output_file, index=False)
                messagebox.showinfo("Success", "Data processed and saved successfully.")
            except Exception as e:
                messagebox.showerror("Error", str(e))
        else:
            messagebox.showerror("Error", "Please check the file paths.")

    def calculate_mz_range(self, formula: str, charge: int = 1, tolerance_ppm: float = 5.0) -> Tuple[float, float]:
        mz = mass.calculate_mass(formula=formula, charge=charge)
        mz_tolerance = mz * tolerance_ppm / 1e6
        return mz - mz_tolerance, mz + mz_tolerance

    def identify_compound(self, row, charge: int, tolerance_ppm: float):
        natural_D = {
            "resiniferonol": ["C20H28O6", "C20H31O6N1"],
            "daphneresiniferin A": ["C24H30O8", "C24H33O8N1"],
            "daphnegen B": ["C24H30O10", "C24H33O10N1"],
            "daphnediterp A": ["C25H36O9", "C25H39O9N1"],
            "daphnetoxin": ["C27H30O8", "C27H33O8N1"],
            "daphne factor F4": ["C27H32O8", "C27H35O8N1"],
            "orthobenzoate 2": ["C27H34O8", "C27H37O8N1"],
            "12-OH-daphnetoxin": ["C27H30O9", "C27H33O9N1"],
            "daphnegiraldigin": ["C27H32O9", "C27H35O9N1"],
            "1,2α-dihydro-5β-hydroxy-6α,7α-epoxy-resiniferonol-14-benzoate": ["C27H34O9", "C27H37O9N1"],
            "1,2β-dihydro-5β-hydroxy-6α,7α-epoxy-resiniferonol-14-benzoate": ["C27H34O9", "C27H37O9N1"],
            "genkwanine O": ["C27H36O9", "C27H39O9N1"],
            "neogenkwanine A": ["C27H36O9", "C27H39O9N1"],
            "neogenkwanine B": ["C27H36O9", "C27H39O9N1"],
            "neogenkwanine H": ["C27H36O9", "C27H39O9N1"],
            "genkwanine A": ["C27H36O9", "C27H39O9N1"],
            "neogenkwanine I": ["C27H36O9", "C27H39O9N1"],
            "daphneresiniferin B": ["C29H32O8", "C29H35O8N1"],
            "genkwadane F": ["C27H32O10", "C27H35O10N1"],
            "genkwanine I": ["C27H38O10", "C27H41O10N1"],
            "yuahuagine": ["C30H38O8", "C30H41O8N1"],
            "excoecaria factor O1": ["C30H38O8", "C30H41O8N1"],
            "yuanhuakine E": ["C29H32O10", "C29H35O10N1"],
            "yuanhuafine": ["C29H32O10", "C29H35O10N1"],
            "yuanhuaoate F": ["C29H32O10", "C29H35O10N1"],
            "yuanhuapine": ["C29H34O10", "C29H37O10N1"],
            "vesiculosin": ["C30H42O9", "C30H45O9N1"],
            "isovesiculosin": ["C30H42O9", "C30H45O9N1"],
            "yuanhuakine C": ["C30H46O9", "C30H49O9N1"],
            "yuanhuaoate A": ["C30H34O10", "C30H37O10N1"],
            "yuanhuakine F": ["C30H36O10", "C30H39O10N1"],
            "wikstrotoxin B": ["C32H44O8", "C32H47O8N1"],
            "genkwanine L": ["C29H36O11", "C29H39O11N1"],
            "daphgenkin F": ["C31H36O10", "C31H39O10N1"],
            "kirkinine D": ["C32H40O10", "C32H43O10N1"],
            "yuanhuadine": ["C32H42O10", "C32H45O10N1"],
            "tianchaterpene D": ["C32H42O10", "C32H45O10N1"],
            "tianchaterpene E": ["C30H40O9", "C30H43O9N1"],
            "isoyuanhuadine": ["C32H42O10", "C32H45O10N1"],
            "yuanhuamine A": ["C32H42O10", "C32H45O10N1"],
            "5β-hydroxyresiniferonol-6α,7α-epoxy-12β-acetoxy-9,13,14-ortho-2E-decenoate": ["C32H44O10", "C32H47O10N1"],
            "yuanhuaoate C": ["C31H42O11", "C31H45O11N1"],
            "genkwanine M": ["C34H38O9", "C34H41O9N1"],
            "genkwanine N": ["C34H38O9", "C34H41O9N1"],
            "yuanhuahine": ["C33H44O10", "C33H47O10N1"],
            "genkwadaphnin": ["C34H34O10", "C34H37O10N1"],
            "tianchaterpene A": ["C32H44O11", "C32H47O11N1"],
            "yuanhuatine": ["C34H36O10", "C34H39O10N1"],
            "genkwadane I": ["C32H44O11", "C32H47O11N1"],
            "genkwanine VIII": ["C34H40O10", "C34H43O10N1"],
            "yuanhuakine D": ["C34H40O10", "C34H43O10N1"],
            "genkwanine D": ["C34H40O10", "C34H43O10N1"],
            "5β-Hydroxyresiniferonol-6α,7α-epoxy- 9,13,14-ortho-2E-hexadecenoate": ["C36H54O8", "C36H57O8N1"],
            "genkwadane D": ["C34H46O10", "C34H49O10N1"],
            "yuahualine": ["C34H46O10", "C34H49O10N1"],
            "stelleralide I": ["C35H50O9", "C35H53O9N1"],
            "acutilonine G": ["C36H40O9", "C36H43O9N1"],
            "wikstroelide N": ["C35H52O9", "C35H55O9N1"],
            "14'-ethyltetrahydrohuratoxin": ["C36H56O8", "C36H59O8N1"],
            "tanguticanine A": ["C35H40O10", "C35H43O10N1"],
            "stelleralide J": ["C34H52O10", "C34H55O10N1"],
            "yuanhuaoate B": ["C34H36O11", "C34H39O11N1"],
            "genkwadane A": ["C34H36O11", "C34H39O11N1"],
            "genkwanine K": ["C34H42O11", "C34H45O11N1"],
            "yuanhuakine G": ["C35H46O10", "C35H49O10N1"],
            "gnidicin": ["C36H36O10", "C36H39O10N1"],
            "yuanhuamine B": ["C35H48O10", "C35H51O10N1"],
            "altadaphnan A": ["C36H38O10", "C36H41O10N1"],
            "acutilonine F": ["C37H46O9", "C37H49O9N1"],
            "wikstroemia factor M1/M2": ["C37H48O9", "C37H51O9N1"],
            "wikstroelide L": ["C36H50O10", "C36H53O10N1"],
            "acetoxyhuratoxin": ["C36H50O10", "C36H53O10N1"],
            "yuanhuamine C": ["C36H50O10", "C36H53O10N1"],
            "altadaphnan C": ["C37H42O10", "C37H45O10N1"],
            "odoratrin": ["C37H42O10", "C37H45O10N1"],
            "daphne factors P1": ["C37H42O10", "C37H45O10N1"],
            "altadaphnan B": ["C37H44O10", "C37H47O10N1"],
            "isoyuanhuacine": ["C37H44O10", "C37H47O10N1"],
            "gnididin": ["C37H44O10", "C37H47O10N1"],
            "hirsein A": ["C37H44O10", "C37H47O10N1"],
            "daphne factors P2": ["C37H44O10", "C37H47O10N1"],
            "tanguticanine C": ["C35H40O12", "C35H43O12N1"],
            "genkwanine C": ["C37H48O10", "C37H51O10N1"],
            "genkwanine E": ["C37H48O10", "C37H51O10N1"],
            "mezerein": ["C38H38O10", "C38H41O10N1"],
            "genkwanine F": ["C37H50O10", "C37H53O10N1"],
            "yuanhuakine B": ["C37H50O10", "C37H53O10N1"],
            "neogenkwanine C": ["C37H50O10", "C37H53O10N1"],
            "neogenkwanine D": ["C37H50O10", "C37H53O10N1"],
            "neogenkwanine E": ["C37H50O10", "C37H53O10N1"],
            "neogenkwanine F": ["C37H50O10", "C37H53O10N1"],
            "genkwanine B": ["C37H50O10", "C37H53O10N1"],
            "genkwanine G": ["C37H52O10", "C37H55O10N1"],
            "wikstroelide B": ["C37H52O10", "C37H55O10N1"],
            "wikstroelide J": ["C36H52O11", "C36H55O11N1"],
            "acutilobin B": ["C37H42O11", "C37H45O11N1"],
            "acutilobin A": ["C37H42O11", "C37H45O11N1"],
            "tanguticanine B": ["C37H44O11", "C37H47O11N1"],
            "genkwadane H": ["C37H46O11", "C37H49O11N1"],
            "tianchaterpene B": ["C37H46O11", "C37H49O11N1"],
            "yuanhuaoate E": ["C37H46O11", "C37H49O11N1"],
            "yuanhuakine A": ["C37H46O11", "C37H49O11N1"],
            "daphgenkin A": ["C37H46O11", "C37H49O11N1"],
            "daphnane-type diterpene ester-7": ["C37H46O11", "C37H49O11N1"],
            "genkwanine J": ["C37H52O11", "C37H55O11N1"],
            "kirkinine": ["C38H56O10", "C38H59O10N1"],
            "acutilobin E": ["C37H38O12", "C37H41O12N1"],
            "tanguticanine D": ["C37H44O12", "C37H47O12N1"],
            "wikstrotoxin C": ["C38H52O11", "C38H55O11N1"],
            "genkwadane G": ["C37H48O12", "C37H51O12N1"],
            "tanguticanine E": ["C38H46O12", "C38H49O12N1"],
            "hirsein B": ["C39H50O11", "C39H53O11N1"],
            "daphnegiraldifin": ["C43H60O9", "C43H63O9N1"],
            "1,2α-dihydro-20-palmitoyldaphnetoxin": ["C43H62O9", "C43H65O9N1"],
            "genkwanine N 20-palmitate": ["C50H68O10", "C50H71O10N1"],
            "genkwadaphnin-20-palmitate": ["C50H64O11", "C50H67O11N1"],
            "wikstroelide C": ["C51H76O11", "C51H79O11N1"],
            "gnidicin-20-palmitate": ["C52H66O11", "C52H69O11N1"],
            "wikstroelide D": ["C52H80O11", "C52H83O11N1"],
            "tanguticacine": ["C53H72O11", "C53H75O11N1"],
            "gnidilatidin 20-palmitate": ["C53H74O11", "C53H77O11N1"],
            "wikstroelide I": ["C53H82O11", "C53H85O11N1"],
            "stelleramacrin B": ["CC34H48O9", "C34H51O9N1"],
            "simplexin (WSC-n)": ["C30H44O8", "C30H47O8N1"],
            "excoecariatoxin (WLI-k)": ["C30H40O8", "C30H43O8N1"],
            "huratoxin (NSC-2)": ["C34H48O8", "C34H51O8N1"],
            "wikstrotoxin A (NSC-3)": ["C35H50O8", "C35H53O8N1"],
            "wikstroelide H (NSC-13)": ["C34H46O10", "C34H49O10N1"],
            "wikstroelide A (NSC-1)": ["C36H50O10", "C36H53O10N1"],
            "yuanhuacine (DOL-C)": ["C37H44O10", "C37H47O10N1"],
            "yuanhuajine (DOL-A)": ["C37H42O10", "C37H45O10N1"],
            "12-O-(E)-cinnamoyl-9,13,14-ortho-(2E,4E)-decadienylidyne-5β,12β-dihydroxyresiniferonol-6α,7α-oxide (DOSOM)": [
                "C39H46O10", "C39H49O10N1"],
            "12-O-(E)-cinnamoyl-9,13,14-ortho-(2E,4E,6E)-decatrienylidyne-5β,12β-dihydroxyresiniferonol-6α,7α-oxide (DOLB)": [
                "C39H44O10", "C39H47O10N1"],
            "daphneodorin D (DOL-H)": ["C39H46O11", "C39H49O11N1"],
            "daphneodorin E (DOL-J)": ["C39H44O11", "C39H47O11N1"],
            "acutilobin D (DOL-I)": ["C40H48O12", "C40H51O12N1"],
            "acutilobin C (DOL-K)": ["C40H46O12", "C40H49O12N1"],
            "gniditrin (DOSO-N)": ["C37H42O10", "C37H45O10N1"],
            "wikstroelide M (NSC-8)": ["C34H50O9", "C34H53O9N1"],
            "stelleralide K (NSC-5)": ["C34H48O11", "C34H51O11N1"],
            "daphneodorin G (DOSO-Q)": ["C39H48O11", "C39H51O11N1"],
            "daphneodorin F (DOL-U)": ["C40H48O13", "C40H51O13N1"],
            "neogenkwanine G (Wlichi-7)": ["C34H40O10", "C34H43O10N1"],
            "daphneodorin H (DOSO-H)": ["C37H48O10", "C37H51O10N1"],
            "genkwanine H (DPSft-2)": ["C34H40O10", "C34H43O10N1"]
        }

        natural_MD = {
            "pimelea factor S6 (TDPE-5)": ["C30H44O8", "C30H47O8N1"],
            "pimelea factor S7 (TDPE-6)": ["C30H44O8", "C30H47O8N1"],
            "wikstroelide F (TDPE-20)": ["C37H48O10", "C37H51O10N1"],
            "daphnepedunin D (TDPE-11)": ["C37H48O11", "C37H51O11N1"],
            "daphnepedunin B (TDPE-10)": ["C30H44O9", "C30H47O9N1"],
            "daphnepedunin A (TDPE-17)": ["C30H44O9", "C30H47O9N1"],
            "daphnepedunin F (TDPE-41)": ["C30H44O9", "C30H47O9N1"],
            "daphnepedunin G (TDPE-42)": ["C30H44O9", "C30H47O9N1"],
            "daphnepedunin C (TDPE-32)": ["C29H42O8", "C29H45O8N1"],
            "stelleralide H (WLI-j)": ["C44H54O12", "C44H57O12N1"],
            "gnidimacrin (DOSO-B)": ["C44H54O12", "C44H57O12N1"],
            "wikstromacrin (WLI-b)": ["C37H50O9", "C37H53O9N1"],
            "pimela factor P2 (WLI-d2)": ["C37H50O9", "C37H53O9N1"],
            "stelleralide G (WIS-12)": ["C51H58O14", "C51H61O14N1"],
            "daphneodorin A (DOSO-K)": ["C53H60O16", "C53H63O16N1"],
            "daphneodorin B (DOSO-E)": ["C55H62O18", "C55H65O18N1"],
            "daphneodorin C (DOL-E2)": ["C55H62O18", "C55H65O18N1"],
            "pimelotide A (NSC-7)": ["C30H42O9", "C30H45O9N1"],
            "pimelotide C (TDPE-4)": ["C30H42O9", "C30H45O9N1"],
            "stelleralide C (TDPE-29)": ["C37H46O11", "C37H49O11N1"],
            "daphnepedunin E (TDPE-28)": ["C37H46O11", "C37H49O11N1"],
            "daphnepedunin F (TDPE-39)": ["C44H50O13", "C44H53O13N1"],
            "daphnopsis factor R1/gnilatimacrin": ["C37H50O9", "C37H53O9N1"],
            "pimelea factor P3": ["C37H50O9", "C37H53O9N1"],
            "genkwadane B": ["C37H50O9", "C37H53O9N1"],
            "linimacrin a": ["C39H52O11", "C39H55O11N1"],
            "linimacrin d": ["C39H52O11", "C39H55O11N1"],
            "stelleralide B": ["C44H54O11", "C44H57O11N1"],
            "kraussianin": ["C37H50O10", "C37H53O10N1"],
            "linimacrin c": ["C37H50O10", "C37H53O10N1"],
            "wikstroelide X": ["C37H50O10", "C37H53O10N1"],
            "stelleramacrin": ["C37H50O11", "C37H53O11N1"],
            "stelleralide A": ["C39H52O12", "C39H55O12N1"],
            "simpleximacrin": ["C46H56O14", "C46H59O14N1"],
            "gnidimacrin-20-palmitate": ["C60H84O13", "C60H87O13N1"],
            "stelleralide F": ["C46H56O14", "C46H59O14N1"],
            "daphnopsis factor R6": ["C37H50O8", "C37H53O8N1"],
            "daphnopsis factor R7": ["C30H46O7", "C30H49O7N1"],
            "wikstroelide E": ["C30H44O8", "C30H47O8N1"],
            "wikstroelide O": ["C37H48O10", "C37H51O10N1"],
            "wikstroelide G": ["C53H78O11", "C53H81O11N1"],
            "wikstroelide K": ["C53H78O11", "C53H81O11N1"],
            "stellejasmin B": ["C37H48O11", "C37H51O11N1"],
            "wikstroelide S": ["C37H52O10", "C37H55O10N1"],
            "genkwadane C": ["C37H52O10", "C37H55O10N1"],
            "wikstroelide T": ["C39H54O10", "C39H57O10N1"],
            "wikstroelide R": ["C37H52O11", "C37H55O11N1"],
            "pimelotide B": ["C32H44O11", "C32H47O11N1"],
            "pimelotide D": ["C32H44O11", "C32H47O11N1"],
            "stelleralide D": ["C53H76O12", "C53H79O12N1"],
            "stelleralide E": ["C55H76O12", "C55H79O12N1"],
            "edgeworthianin D (ECH-1)": ["C36H52O10", "C36H55O10N1"],
            "edgeworthianin E (ECH-2)": ["C40H58O12", "C40H61O12N1"],
            "edgeworthianin A (ECH-3)": ["C36H50O11", "C36H53O11N1"],
            "edgeworthianin B (ECH-4)": ["C40H56O13", "C40H59O13N1"],
            "edgeworthianin C (ECH-5)": ["C43H54O13", "C43H57O13N1"],
            "edgeworthianin G (ECS-6)": ["C41H58O13", "C41H61O13N1"],
            "edgeworthianin F (ECF-7)": ["C38H54O12", "C38H57O12N1"],
            "kirkinine B": ["C36H54O8", "C36H57O8N1"],
            "kirkinine E": ["C36H54O9", "C36H57O9N1"],
            "kirkinine C": ["C38H56O10", "C38H59O10N1"],
            "20-O-hexadecanoylkirkinine B": ["C52H84O9", "C52H87O9N1"],
            "synaptolepis factor K1": ["C36H54O8", "C36H57O8N1"],
        }

        known_compounds = natural_D if row['type'] == 'D' else natural_MD if row['type'] == 'MD' else {}
        identifications = []
        for compound, formulas in known_compounds.items():
            for formula in formulas:
                mz_min, mz_max = self.calculate_mz_range(formula, charge, tolerance_ppm)
                if mz_min <= row['row m/z'] <= mz_max:
                    identifications.append(compound)
        if identifications:
            return 'possibly known', ', '.join(set(identifications))
        else:
            return 'possibly undescribed', 'none'

    def calculate_formula_difference(self, formula1: str, formula2: str) -> dict:
        comp1 = mass.Composition(formula1)
        comp2 = mass.Composition(formula2)
        mw1 = mass.calculate_mass(composition=comp1)
        mw2 = mass.calculate_mass(composition=comp2)
        if mw1 > mw2:
            diff = {elem: comp1.get(elem, 0) - comp2.get(elem, 0) for elem in ('C', 'H', 'O', 'N')}
        else:
            diff = {elem: comp2.get(elem, 0) - comp1.get(elem, 0) for elem in ('C', 'H', 'O', 'N')}
        return diff

    def h2o_or_nh3(self, diff):
        if diff == {'C': 0, 'H': 2, 'O': 1, 'N': 0}:
            return "H2O"
        elif diff == {'C': 0, 'H': 3, 'O': 0, 'N': 1}:
            return "NH3"
        return None

    def create_rt_alignment(self, df):
        rt_alignment_dict = {}
        duplicated_rts = df['row retention time'].duplicated(keep=False)
        for rt in df[duplicated_rts]['row retention time'].unique():
            indices = df.index[df['row retention time'] == rt].tolist()
            row_ids = df.loc[indices, 'row ID'].astype(str).tolist()
            for index in indices:
                rt_alignment_dict[index] = ', '.join(row_ids)
        return [rt_alignment_dict.get(index, 'none') for index in df.index]

    def add_annotations(self, df):
        df['adduct ion'] = df['elemental composition'].apply(lambda x: '[M+NH4]+' if 'N' in x else '[M+H]+')
        duplicated_rts = df['row retention time'].duplicated(keep=False)
        for rt in df[duplicated_rts]['row retention time'].unique():
            indices = df.index[df['row retention time'] == rt].tolist()
            for i in range(len(indices) - 1):
                for j in range(i + 1, len(indices)):
                    diff = self.calculate_formula_difference(df.at[indices[i], 'elemental composition'],
                                                            df.at[indices[j], 'elemental composition'])
                    if self.h2o_or_nh3(diff) == "H2O":
                        if (mass.Composition(df.at[indices[i], 'elemental composition']).get('O', 0)
                                < mass.Composition(df.at[indices[j], 'elemental composition']).get('O', 0)):
                            df.at[indices[i], 'adduct ion'] = '[M+H-H2O]+'
                        else:
                            df.at[indices[j], 'adduct ion'] = '[M+H-H2O]+'
        return df

    def compare_peaks_and_mark(self, df):
        same_peak_dict = {}
        c_part_dict = {}
        grouped = df.groupby('row retention time')
        for name, group in grouped:
            if len(group) > 1:
                ids = group.index
                for i in range(len(ids)):
                    for j in range(i + 1, len(ids)):
                        type1 = group.at[ids[i], 'type']
                        type2 = group.at[ids[j], 'type']
                        formula1 = group.at[ids[i], 'elemental composition']
                        formula2 = group.at[ids[j], 'elemental composition']
                        diff = self.calculate_formula_difference(formula1, formula2)
                        result = self.h2o_or_nh3(diff)
                        if result in ["H2O", "NH3"]:
                            if ids[i] in same_peak_dict:
                                same_peak_dict[ids[i]].add(group.at[ids[j], 'row ID'])
                            else:
                                same_peak_dict[ids[i]] = {group.at[ids[i], 'row ID'], group.at[ids[j], 'row ID']}
                            if ids[j] in same_peak_dict:
                                same_peak_dict[ids[j]].add(group.at[ids[i], 'row ID'])
                            else:
                                same_peak_dict[ids[j]] = {group.at[ids[i], 'row ID'], group.at[ids[j], 'row ID']}
                        if result == "H2O" and type1 == "D" and type2 == "D":
                            c_part_dict[ids[i]] = "9, 13, 14-triol"
                            c_part_dict[ids[j]] = "9, 13, 14-triol"
        for key, value in same_peak_dict.items():
            same_peak_dict[key] = ', '.join(sorted(str(v) for v in value))
        for idx in df.index:
            if idx not in same_peak_dict:
                same_peak_dict[idx] = 'none'
        for idx in df.index:
            if idx not in c_part_dict:
                c_part_dict[idx] = "9, 13, 14-orthoester"
        df['same peak'] = pd.Series(same_peak_dict)
        df['C part'] = pd.Series(c_part_dict)
        return df

    def process_dataframe(self, df: pd.DataFrame, charge: int, tolerance_ppm: float) -> pd.DataFrame:
        df['identification'], df['same composition'] = zip(*df.apply(self.identify_compound, args=(charge, tolerance_ppm), axis=1))
        df = self.add_annotations(df)
        ec_index = df.columns.get_loc('elemental composition') + 1
        df.insert(ec_index, 'adduct ion', df.pop('adduct ion'))
        df = self.compare_peaks_and_mark(df)
        ec_index = df.columns.get_loc('row retention time') + 1
        df.insert(ec_index, 'same peak', df.pop('same peak'))
        df['row retention time'] = df['row retention time'].round(2)
        peak_area_column = next(col for col in df.columns if 'Peak area' in col)
        df.sort_values(by=peak_area_column, ascending=False, inplace=True)
        return df


    '''neutral_loss_extractor_tab'''
    def create_neutral_loss_extractor_tab(self):
        neutral_loss_extractor_tab = self.tab("4.Neutral Loss Extract")

        # Create a canvas
        canvas = tk.Canvas(neutral_loss_extractor_tab)
        canvas.pack(side="left", fill="both", expand=True)

        # Add a scrollbar
        scrollbar = tk.Scrollbar(neutral_loss_extractor_tab, orient="vertical", command=canvas.yview)
        scrollbar.pack(side="right", fill="y")

        # Configure canvas with scrollbar
        canvas.configure(yscrollcommand=scrollbar.set)
        canvas.bind("<Configure>", lambda e: canvas.configure(scrollregion=canvas.bbox("all")))

        # Create a frame inside the canvas
        frame_inside_canvas = ctk.CTkFrame(canvas)
        canvas.create_window((0, 0), window=frame_inside_canvas, anchor='nw')

        # Add the existing widgets to the frame inside the canvas
        # CTkFrame of Import MGF
        frame_ImportMGF = ctk.CTkFrame(frame_inside_canvas, fg_color="transparent")
        frame_ImportMGF.pack(fill='x', padx=10, pady=5)
        label_NLextractor_mgf = ctk.CTkLabel(frame_ImportMGF, text="Import MGF file:", font=('Arial', 14, 'bold'))
        label_NLextractor_mgf.pack(side='left', padx=(0, 10))
        self.entry_NLextractor_mgf_path = ctk.CTkEntry(frame_ImportMGF, placeholder_text="Input file path", width=548)
        self.entry_NLextractor_mgf_path.pack(side='left', padx=(0, 5))
        self.button_NLextractor_browse = ctk.CTkButton(frame_ImportMGF, text="Browse file", command=self.browse_NLextractor_mgf, fg_color="#6598D3", hover_color="#5578B3")
        self.button_NLextractor_browse.pack(padx=10, pady=5, anchor='e')

        # CTkFrame of Import CSV
        frame_ImportCSV = ctk.CTkFrame(frame_inside_canvas, fg_color="transparent")
        frame_ImportCSV.pack(fill='x', padx=10, pady=5)
        label_NLextractor_csv = ctk.CTkLabel(frame_ImportCSV, text="Import CSV file:", font=('Arial', 14, 'bold'))
        label_NLextractor_csv.pack(side='left', padx=(0, 10))
        self.entry_NLextractor_csv_path = ctk.CTkEntry(frame_ImportCSV, placeholder_text="Input CSV file path", width=549)
        self.entry_NLextractor_csv_path.pack(side='left', padx=(0, 5))
        self.button_NLcsv_browse = ctk.CTkButton(frame_ImportCSV, text="Browse file", command=self.browse_NLextractor_csv, fg_color="#6598D3", hover_color="#5578B3")
        self.button_NLcsv_browse.pack(padx=10, pady=5, anchor='e')

        # CTkFrame of Import DB
        frame_ImportDB = ctk.CTkFrame(frame_inside_canvas, fg_color="transparent")
        frame_ImportDB.pack(fill='x', padx=10, pady=5)
        label_db = ctk.CTkLabel(frame_ImportDB, text="Import formula database:", font=('Arial', 14, 'bold'))
        label_db.pack(side='left', padx=(0, 10))
        self.entry_NLdb_path = ctk.CTkEntry(frame_ImportDB, placeholder_text="Input database file path", width=485)
        self.entry_NLdb_path.pack(side='left', padx=(0, 5))
        self.button_NLdb_browse = ctk.CTkButton(frame_ImportDB, text="Browse file", command=self.browse_NLdb, fg_color="#6598D3", hover_color="#5578B3")
        self.button_NLdb_browse.pack(padx=10, pady=5, anchor='e')

        # CTkFrame of Mass Tolerance (ppm) and Charge
        frame_charge_ppm = ctk.CTkFrame(frame_inside_canvas, fg_color="transparent")
        frame_charge_ppm.pack(fill='x', padx=10, pady=5)
        # Charge
        label_charge = ctk.CTkLabel(frame_charge_ppm, text="・Charge:", font=('Arial', 14))
        label_charge.pack(side='left', padx=(0, 10))
        self.entry_charge = ctk.CTkEntry(frame_charge_ppm, placeholder_text="Input charge (e.g., +1,-1)", width=278)
        self.entry_charge.pack(side='left', padx=(0, 20))
        # Mass Tolerance (ppm)
        label_ppm = ctk.CTkLabel(frame_charge_ppm, text="・Mass Tolerance (ppm):", font=('Arial', 14))
        label_ppm.pack(side='left', padx=(0, 10))
        self.entry_ppm = ctk.CTkEntry(frame_charge_ppm, placeholder_text="Input ppm (e.g., 5)", width=278)
        self.entry_ppm.pack(side='left')

        # Targerloss 1
        frame_targetloss1 = ctk.CTkFrame(frame_inside_canvas, fg_color="transparent")
        frame_targetloss1.pack(fill='x', padx=10, pady=5)
        label_targetloss1 = ctk.CTkLabel(frame_targetloss1, text="  Target Loss 1 (or set):", font=('Arial', 14, 'bold'))
        label_targetloss1.pack(side='left', padx=(0, 10))
        self.entry_targetloss1 = ctk.CTkEntry(frame_targetloss1, placeholder_text="Loss 1 (e.g., C3H4O2)", width=345)
        self.entry_targetloss1.pack(side='left', padx=(0, 10))
        # Loss name 1
        label_lossname1 = ctk.CTkLabel(frame_targetloss1, text="Name 1:", font=('Arial', 14, 'bold'))
        label_lossname1.pack(side='left', padx=(0, 10))
        self.entry_lossname1 = ctk.CTkEntry(frame_targetloss1, placeholder_text="e.g., C3 loss", width=150)
        self.entry_lossname1.pack(side='left', padx=(0, 10))

        # Targerloss 1_1
        frame_targetloss1_1 = ctk.CTkFrame(frame_inside_canvas, fg_color="transparent")
        frame_targetloss1_1.pack(fill='x', padx=10, pady=5)
        # Carbon Range
        label_carbonRange1_1_max = ctk.CTkLabel(frame_targetloss1_1, text="・Carbon Range 1  From:", font=('Arial', 14))
        label_carbonRange1_1_max.pack(side='left', padx=(0, 10))
        self.entry_carbonRange1_1_max = ctk.CTkEntry(frame_targetloss1_1, placeholder_text="C (max), e.g., C20", width=150)
        self.entry_carbonRange1_1_max.pack(side='left', padx=(0, 10))
        # Carbon Range
        label_carbonRange1_1_min = ctk.CTkLabel(frame_targetloss1_1, text="To:", font=('Arial', 14))
        label_carbonRange1_1_min.pack(side='left', padx=(0, 10))
        self.entry_carbonRange1_1_min = ctk.CTkEntry(frame_targetloss1_1, placeholder_text="C (min), e.g., C17", width=150)
        self.entry_carbonRange1_1_min.pack(side='left', padx=(0, 10))

        # Targerloss 1_2
        frame_targetloss1_2 = ctk.CTkFrame(frame_inside_canvas, fg_color="transparent")
        frame_targetloss1_2.pack(fill='x', padx=10, pady=5)
        # Carbon Range
        label_carbonRange1_2_max = ctk.CTkLabel(frame_targetloss1_2, text="・Carbon Range 2  From:", font=('Arial', 14))
        label_carbonRange1_2_max.pack(side='left', padx=(0, 10))
        self.entry_carbonRange1_2_max = ctk.CTkEntry(frame_targetloss1_2, placeholder_text="C (max), e.g., C30", width=150)
        self.entry_carbonRange1_2_max.pack(side='left', padx=(0, 10))
        # Carbon Range
        label_carbonRange1_2_min = ctk.CTkLabel(frame_targetloss1_2, text="To:", font=('Arial', 14))
        label_carbonRange1_2_min.pack(side='left', padx=(0, 10))
        self.entry_carbonRange1_2_min = ctk.CTkEntry(frame_targetloss1_2, placeholder_text="C (min), e.g., C27", width=150)
        self.entry_carbonRange1_2_min.pack(side='left', padx=(0, 10))

        # Targerloss 2
        frame_targetloss2 = ctk.CTkFrame(frame_inside_canvas, fg_color="transparent")
        frame_targetloss2.pack(fill='x', padx=10, pady=5)
        label_targetloss2 = ctk.CTkLabel(frame_targetloss2, text="  Target Loss 2 (or set):", font=('Arial', 14, 'bold'))
        label_targetloss2.pack(side='left', padx=(0, 10))
        self.entry_targetloss2 = ctk.CTkEntry(frame_targetloss2, placeholder_text="Loss 2 (e.g., C10H18O2, C10H16O2)", width=345)
        # C10H18O2, C10H16O2, C10H14O2, C10H12O2, C10H10O2
        self.entry_targetloss2.pack(side='left', padx=(0, 10))
        # Loss name 2
        label_lossname2 = ctk.CTkLabel(frame_targetloss2, text="Name 2:", font=('Arial', 14, 'bold'))
        label_lossname2.pack(side='left', padx=(0, 10))
        self.entry_lossname2 = ctk.CTkEntry(frame_targetloss2, placeholder_text="e.g., C10 loss", width=150)
        self.entry_lossname2.pack(side='left', padx=(0, 10))

        # Targerloss 2_1
        frame_targetloss2_1 = ctk.CTkFrame(frame_inside_canvas, fg_color="transparent")
        frame_targetloss2_1.pack(fill='x', padx=10, pady=5)
        # Carbon Range
        label_carbonRange2_1_max = ctk.CTkLabel(frame_targetloss2_1, text="・Carbon Range 1  From:", font=('Arial', 14))
        label_carbonRange2_1_max.pack(side='left', padx=(0, 10))
        self.entry_carbonRange2_1_max = ctk.CTkEntry(frame_targetloss2_1, placeholder_text="C (max), e.g., C30", width=150)
        self.entry_carbonRange2_1_max.pack(side='left', padx=(0, 10))
        # Carbon Range
        label_carbonRange2_1_min = ctk.CTkLabel(frame_targetloss2_1, text="To:", font=('Arial', 14))
        label_carbonRange2_1_min.pack(side='left', padx=(0, 10))
        self.entry_carbonRange2_1_min = ctk.CTkEntry(frame_targetloss2_1, placeholder_text="C (min), e.g., C20", width=150)
        self.entry_carbonRange2_1_min.pack(side='left', padx=(0, 10))

        # Targerloss 2_2
        frame_targetloss2_2 = ctk.CTkFrame(frame_inside_canvas, fg_color="transparent")
        frame_targetloss2_2.pack(fill='x', padx=10, pady=5)
        # Carbon Range
        label_carbonRange2_2_max = ctk.CTkLabel(frame_targetloss2_2, text="・Carbon Range 2  From:", font=('Arial', 14))
        label_carbonRange2_2_max.pack(side='left', padx=(0, 10))
        self.entry_carbonRange2_2_max = ctk.CTkEntry(frame_targetloss2_2, placeholder_text="C (max)", width=150)
        self.entry_carbonRange2_2_max.pack(side='left', padx=(0, 10))
        # Carbon Range
        label_carbonRange2_2_min = ctk.CTkLabel(frame_targetloss2_2, text="To:", font=('Arial', 14))
        label_carbonRange2_2_min.pack(side='left', padx=(0, 10))
        self.entry_carbonRange2_2_min = ctk.CTkEntry(frame_targetloss2_2, placeholder_text="C (min)", width=150)
        self.entry_carbonRange2_2_min.pack(side='left', padx=(0, 10))

        # Targerloss 3
        frame_targetloss3 = ctk.CTkFrame(frame_inside_canvas, fg_color="transparent")
        frame_targetloss3.pack(fill='x', padx=10, pady=5)
        label_targetloss3 = ctk.CTkLabel(frame_targetloss3, text="  Target Loss 3 (or set):", font=('Arial', 14, 'bold'))
        label_targetloss3.pack(side='left', padx=(0, 10))
        self.entry_targetloss3 = ctk.CTkEntry(frame_targetloss3, placeholder_text="Loss 3 (e.g., C1H0O2)", width=345)
        self.entry_targetloss3.pack(side='left', padx=(0, 10))
        # Loss name 2
        label_lossname3 = ctk.CTkLabel(frame_targetloss3, text="Name 3:", font=('Arial', 14, 'bold'))
        label_lossname3.pack(side='left', padx=(0, 10))
        self.entry_lossname3 = ctk.CTkEntry(frame_targetloss3, placeholder_text="e.g., CO2 loss", width=150)
        self.entry_lossname3.pack(side='left', padx=(0, 10))

        # Targerloss 3_1
        frame_targetloss3_1 = ctk.CTkFrame(frame_inside_canvas, fg_color="transparent")
        frame_targetloss3_1.pack(fill='x', padx=10, pady=5)
        # Carbon Range
        label_carbonRange3_1_max = ctk.CTkLabel(frame_targetloss3_1, text="・Carbon Range 1  From:", font=('Arial', 14))
        label_carbonRange3_1_max.pack(side='left', padx=(0, 10))
        self.entry_carbonRange3_1_max = ctk.CTkEntry(frame_targetloss3_1, placeholder_text="C (max), e.g., C30", width=150)
        self.entry_carbonRange3_1_max.pack(side='left', padx=(0, 10))
        # Carbon Range
        label_carbonRange3_1_min = ctk.CTkLabel(frame_targetloss3_1, text="To:", font=('Arial', 14))
        label_carbonRange3_1_min.pack(side='left', padx=(0, 10))
        self.entry_carbonRange3_1_min = ctk.CTkEntry(frame_targetloss3_1, placeholder_text="C (min), e.g., C29", width=150)
        self.entry_carbonRange3_1_min.pack(side='left', padx=(0, 10))

        # Targerloss 3_2
        frame_targetloss3_2 = ctk.CTkFrame(frame_inside_canvas, fg_color="transparent")
        frame_targetloss3_2.pack(fill='x', padx=10, pady=5)
        # Carbon Range
        label_carbonRange3_2_max = ctk.CTkLabel(frame_targetloss3_2, text="・Carbon Range 2  From:", font=('Arial', 14))
        label_carbonRange3_2_max.pack(side='left', padx=(0, 10))
        self.entry_carbonRange3_2_max = ctk.CTkEntry(frame_targetloss3_2, placeholder_text="C (max)", width=150)
        self.entry_carbonRange3_2_max.pack(side='left', padx=(0, 10))
        # Carbon Range
        label_carbonRange3_2_min = ctk.CTkLabel(frame_targetloss3_2, text="To:", font=('Arial', 14))
        label_carbonRange3_2_min.pack(side='left', padx=(0, 10))
        self.entry_carbonRange3_2_min = ctk.CTkEntry(frame_targetloss3_2, placeholder_text="C (min)", width=150)
        self.entry_carbonRange3_2_min.pack(side='left', padx=(0, 10))

        # Targerloss 4
        frame_targetloss4 = ctk.CTkFrame(frame_inside_canvas, fg_color="transparent")
        frame_targetloss4.pack(fill='x', padx=10, pady=5)
        label_targetloss4 = ctk.CTkLabel(frame_targetloss4, text="  Target Loss 4 (or set):", font=('Arial', 14, 'bold'))
        label_targetloss4.pack(side='left', padx=(0, 10))
        self.entry_targetloss4 = ctk.CTkEntry(frame_targetloss4, placeholder_text="Loss 4 (e.g., C2H4O2)", width=345)
        self.entry_targetloss4.pack(side='left', padx=(0, 10))
        # Loss name 2
        label_lossname4 = ctk.CTkLabel(frame_targetloss4, text="Name 4:", font=('Arial', 14, 'bold'))
        label_lossname4.pack(side='left', padx=(0, 10))
        self.entry_lossname4 = ctk.CTkEntry(frame_targetloss4, placeholder_text="e.g., C2 loss", width=150)
        self.entry_lossname4.pack(side='left', padx=(0, 10))

        # Targerloss 4_1
        frame_targetloss4_1 = ctk.CTkFrame(frame_inside_canvas, fg_color="transparent")
        frame_targetloss4_1.pack(fill='x', padx=10, pady=5)
        # Carbon Range
        label_carbonRange4_1_max = ctk.CTkLabel(frame_targetloss4_1, text="・Carbon Range 1  From:", font=('Arial', 14))
        label_carbonRange4_1_max.pack(side='left', padx=(0, 10))
        self.entry_carbonRange4_1_max = ctk.CTkEntry(frame_targetloss4_1, placeholder_text="C (max), e.g., C32", width=150)
        self.entry_carbonRange4_1_max.pack(side='left', padx=(0, 10))
        # Carbon Range
        label_carbonRange4_1_min = ctk.CTkLabel(frame_targetloss4_1, text="To:", font=('Arial', 14))
        label_carbonRange4_1_min.pack(side='left', padx=(0, 10))
        self.entry_carbonRange4_1_min = ctk.CTkEntry(frame_targetloss4_1, placeholder_text="C (min), e.g., C30", width=150)
        self.entry_carbonRange4_1_min.pack(side='left', padx=(0, 10))

        # Targerloss 4_2
        frame_targetloss4_2 = ctk.CTkFrame(frame_inside_canvas, fg_color="transparent")
        frame_targetloss4_2.pack(fill='x', padx=10, pady=5)
        # Carbon Range
        label_carbonRange4_2_max = ctk.CTkLabel(frame_targetloss4_2, text="・Carbon Range 2  From:", font=('Arial', 14))
        label_carbonRange4_2_max.pack(side='left', padx=(0, 10))
        self.entry_carbonRange4_2_max = ctk.CTkEntry(frame_targetloss4_2, placeholder_text="C (max)e.g., C22", width=150)
        self.entry_carbonRange4_2_max.pack(side='left', padx=(0, 10))
        # Carbon Range
        label_carbonRange4_2_min = ctk.CTkLabel(frame_targetloss4_2, text="To:", font=('Arial', 14))
        label_carbonRange4_2_min.pack(side='left', padx=(0, 10))
        self.entry_carbonRange4_2_min = ctk.CTkEntry(frame_targetloss4_2, placeholder_text="C (min)e.g., C20", width=150)
        self.entry_carbonRange4_2_min.pack(side='left', padx=(0, 10))

        # Targerloss 5
        frame_targetloss5 = ctk.CTkFrame(frame_inside_canvas, fg_color="transparent")
        frame_targetloss5.pack(fill='x', padx=10, pady=5)
        label_targetloss5 = ctk.CTkLabel(frame_targetloss5, text="  Target Loss 5 (or set):", font=('Arial', 14, 'bold'))
        label_targetloss5.pack(side='left', padx=(0, 10))
        self.entry_targetloss5 = ctk.CTkEntry(frame_targetloss5, placeholder_text="Loss 5 (e.g., C7H6O2)", width=345)
        self.entry_targetloss5.pack(side='left', padx=(0, 10))
        # Loss name 2
        label_lossname5 = ctk.CTkLabel(frame_targetloss5, text="Name 5:", font=('Arial', 14, 'bold'))
        label_lossname5.pack(side='left', padx=(0, 10))
        self.entry_lossname5 = ctk.CTkEntry(frame_targetloss5, placeholder_text="e.g., C7 loss", width=150)
        self.entry_lossname5.pack(side='left', padx=(0, 10))

        # Targerloss 5_1
        frame_targetloss5_1 = ctk.CTkFrame(frame_inside_canvas, fg_color="transparent")
        frame_targetloss5_1.pack(fill='x', padx=10, pady=5)
        # Carbon Range
        label_carbonRange5_1_max = ctk.CTkLabel(frame_targetloss5_1, text="・Carbon Range 1  From:", font=('Arial', 14))
        label_carbonRange5_1_max.pack(side='left', padx=(0, 10))
        self.entry_carbonRange5_1_max = ctk.CTkEntry(frame_targetloss5_1, placeholder_text="C (max), e.g., C37", width=150)
        self.entry_carbonRange5_1_max.pack(side='left', padx=(0, 10))
        # Carbon Range
        label_carbonRange5_1_min = ctk.CTkLabel(frame_targetloss5_1, text="To:", font=('Arial', 14))
        label_carbonRange5_1_min.pack(side='left', padx=(0, 10))
        self.entry_carbonRange5_1_min = ctk.CTkEntry(frame_targetloss5_1, placeholder_text="C (min), e.g., C30", width=150)
        self.entry_carbonRange5_1_min.pack(side='left', padx=(0, 10))

        # Targerloss 5_2
        frame_targetloss5_2 = ctk.CTkFrame(frame_inside_canvas, fg_color="transparent")
        frame_targetloss5_2.pack(fill='x', padx=10, pady=5)
        # Carbon Range
        label_carbonRange5_2_max = ctk.CTkLabel(frame_targetloss5_2, text="・Carbon Range 2  From:", font=('Arial', 14))
        label_carbonRange5_2_max.pack(side='left', padx=(0, 10))
        self.entry_carbonRange5_2_max = ctk.CTkEntry(frame_targetloss5_2, placeholder_text="C (max)e.g., C27", width=150)
        self.entry_carbonRange5_2_max.pack(side='left', padx=(0, 10))
        # Carbon Range
        label_carbonRange5_2_min = ctk.CTkLabel(frame_targetloss5_2, text="To:", font=('Arial', 14))
        label_carbonRange5_2_min.pack(side='left', padx=(0, 10))
        self.entry_carbonRange5_2_min = ctk.CTkEntry(frame_targetloss5_2, placeholder_text="C (min)e.g., C20", width=150)
        self.entry_carbonRange5_2_min.pack(side='left', padx=(0, 10))

        # Targerloss 6
        frame_targetloss6 = ctk.CTkFrame(frame_inside_canvas, fg_color="transparent")
        frame_targetloss6.pack(fill='x', padx=10, pady=5)
        label_targetloss6 = ctk.CTkLabel(frame_targetloss6, text="  Target Loss 6 (or set):",
                                         font=('Arial', 14, 'bold'))
        label_targetloss6.pack(side='left', padx=(0, 10))
        self.entry_targetloss6 = ctk.CTkEntry(frame_targetloss6, placeholder_text="Loss 6 (e.g., C14H24O2)", width=345)
        self.entry_targetloss6.pack(side='left', padx=(0, 10))
        # Loss name 2
        label_lossname6 = ctk.CTkLabel(frame_targetloss6, text="Name 6:", font=('Arial', 14, 'bold'))
        label_lossname6.pack(side='left', padx=(0, 10))
        self.entry_lossname6 = ctk.CTkEntry(frame_targetloss6, placeholder_text="e.g., C14 loss", width=150)
        self.entry_lossname6.pack(side='left', padx=(0, 10))

        # Targerloss 6_1
        frame_targetloss6_1 = ctk.CTkFrame(frame_inside_canvas, fg_color="transparent")
        frame_targetloss6_1.pack(fill='x', padx=10, pady=5)
        # Carbon Range
        label_carbonRange6_1_max = ctk.CTkLabel(frame_targetloss6_1, text="・Carbon Range 1  From:", font=('Arial', 14))
        label_carbonRange6_1_max.pack(side='left', padx=(0, 10))
        self.entry_carbonRange6_1_max = ctk.CTkEntry(frame_targetloss6_1, placeholder_text="C (max), e.g., C34",
                                                     width=150)
        self.entry_carbonRange6_1_max.pack(side='left', padx=(0, 10))
        # Carbon Range
        label_carbonRange6_1_min = ctk.CTkLabel(frame_targetloss6_1, text="To:", font=('Arial', 14))
        label_carbonRange6_1_min.pack(side='left', padx=(0, 10))
        self.entry_carbonRange6_1_min = ctk.CTkEntry(frame_targetloss6_1, placeholder_text="C (min), e.g., C20",
                                                     width=150)
        self.entry_carbonRange6_1_min.pack(side='left', padx=(0, 10))

        # Targerloss 6_2
        frame_targetloss6_2 = ctk.CTkFrame(frame_inside_canvas, fg_color="transparent")
        frame_targetloss6_2.pack(fill='x', padx=10, pady=5)
        # Carbon Range
        label_carbonRange6_2_max = ctk.CTkLabel(frame_targetloss6_2, text="・Carbon Range 2  From:", font=('Arial', 14))
        label_carbonRange6_2_max.pack(side='left', padx=(0, 10))
        self.entry_carbonRange6_2_max = ctk.CTkEntry(frame_targetloss6_2, placeholder_text="C (max)e.g., C36",
                                                     width=150)
        self.entry_carbonRange6_2_max.pack(side='left', padx=(0, 10))
        # Carbon Range
        label_carbonRange6_2_min = ctk.CTkLabel(frame_targetloss6_2, text="To:", font=('Arial', 14))
        label_carbonRange6_2_min.pack(side='left', padx=(0, 10))
        self.entry_carbonRange6_2_min = ctk.CTkEntry(frame_targetloss6_2, placeholder_text="C (min)e.g., C22",
                                                     width=150)
        self.entry_carbonRange6_2_min.pack(side='left', padx=(0, 10))

        # CTkFrame of Import Output
        frame_Output = ctk.CTkFrame(frame_inside_canvas, fg_color="transparent")
        frame_Output.pack(fill='x', padx=10, pady=5)
        label_output = ctk.CTkLabel(frame_Output, text="Export Result Files:", font=('Arial', 14, 'bold'))
        label_output.pack(side='left', padx=(0, 10))
        self.entry_NLoutput_path = ctk.CTkEntry(frame_Output, placeholder_text="Output csv file path, e.g. C:/Users/xxx/", width=314)
        self.entry_NLoutput_path.pack(side='left', padx=(0, 5))
        self.entry_NLoutput_csv_filename = ctk.CTkEntry(frame_Output, placeholder_text="Output csv file name, e.g. xxx.csv", width=202)
        self.entry_NLoutput_csv_filename.pack(side='left', padx=(0, 5))
        self.button_NLextract = ctk.CTkButton(frame_Output, text="Extract Neutral Loss", command=self.NLextract, fg_color="#6598D3", hover_color="#5578B3")
        self.button_NLextract.pack(padx=10, pady=5, anchor='e')

        # Progress Bar
        self.progress_bar = ctk.CTkProgressBar(frame_Output, width=0)
        self.progress_bar.pack(fill='x', padx=10)
        self.progress_bar.set(0)

    def browse_NLextractor_mgf(self):
        path = filedialog.askopenfilename(filetypes=[("MGF files", "*.mgf")])
        if path:
            self.entry_NLextractor_mgf_path.delete(0, tk.END)
            self.entry_NLextractor_mgf_path.insert(0, path)

    def browse_NLextractor_csv(self):
        path = filedialog.askopenfilename(filetypes=[("CSV files", "*.csv")])
        if path:
            self.entry_NLextractor_csv_path.delete(0, tk.END)
            self.entry_NLextractor_csv_path.insert(0, path)

    def browse_NLdb(self):
        path = filedialog.askopenfilename(filetypes=[("Database files", "*.db")])
        if path:
            self.entry_NLdb_path.delete(0, tk.END)
            self.entry_NLdb_path.insert(0, path)

    def NLextract(self):
        def parse_formula(formula):
            try:
                if formula and formula.lower() != 'none':
                    return Composition(formula)
                else:
                    raise ValueError("Empty or placeholder formula given.")
            except Exception as e:
                print(f"Error parsing formula '{formula}': {str(e)}")
                raise

        try:
            target_loss1 = [parse_formula(x.strip()) for x in self.entry_targetloss1.get().split(',') if x.strip()]
            target_loss2 = [parse_formula(x.strip()) for x in self.entry_targetloss2.get().split(',') if x.strip()]
            target_loss3 = [parse_formula(x.strip()) for x in self.entry_targetloss3.get().split(',') if x.strip()]
            target_loss4 = [parse_formula(x.strip()) for x in self.entry_targetloss4.get().split(',') if x.strip()]
            target_loss5 = [parse_formula(x.strip()) for x in self.entry_targetloss5.get().split(',') if x.strip()]
            target_loss6 = [parse_formula(x.strip()) for x in self.entry_targetloss6.get().split(',') if x.strip()]

            charge = int(self.entry_charge.get().strip())
            ppm = float(self.entry_ppm.get().strip())
            output_dir = self.entry_NLoutput_path.get().strip()
            output_filename = self.entry_NLoutput_csv_filename.get().strip()
            full_output_path = os.path.join(output_dir, output_filename)

            self.process_files(
                self.entry_NLextractor_csv_path.get(),
                self.entry_NLextractor_mgf_path.get(),
                full_output_path,
                self.entry_NLdb_path.get(),
                target_loss1, target_loss2, target_loss3, target_loss4, target_loss5, target_loss6,
                charge, ppm
            )
            messagebox.showinfo("Process Complete", "Neutral losses have been successfully extracted.")
        except ValueError as ve:
            messagebox.showerror("Error", str(ve))
        except Exception as e:
            messagebox.showerror("Error", f"Failed to process data: {str(e)}")

    def mz_range(self, formula: str, charge: int, ppm: float) -> Tuple[float, float]:
        mz = mass.calculate_mass(formula=formula, charge=charge)
        mz_min = mz - mz * ppm / 1e6
        mz_max = mz + mz * ppm / 1e6
        return abs(mz_min), abs(mz_max)

    def read_formulas_from_NLdb(self, db_path: str, charge: int, ppm: float) -> Dict[str, Tuple[float, float]]:
        conn = sqlite3.connect(db_path)
        cursor = conn.cursor()
        cursor.execute("SELECT formula FROM compounds")
        formulas = cursor.fetchall()
        conn.close()
        return {formula[0]: self.mz_range(formula[0], charge, ppm) for formula in formulas}

    def find_matching_NLformulas(self, mz_value: float, formula_ranges: dict) -> str:
        matching_formulas = [
            formula for formula, (mz_min, mz_max) in formula_ranges.items()
            if mz_min <= mz_value <= mz_max
        ]
        print(f"Matching formulas for mz {mz_value}: {matching_formulas}")
        return matching_formulas[0] if len(matching_formulas) == 1 else 'none'

    def specific_formula_difference(self, formula1: str, formula2: str) -> Dict[str, int]:
        comp1 = mass.Composition(formula1)
        comp2 = mass.Composition(formula2)
        differences = {
            'C': abs(comp2.get('C', 0) - comp1.get('C', 0)),
            'H': abs(comp2.get('H', 0) - comp1.get('H', 0)),
            'O': abs(comp2.get('O', 0) - comp1.get('O', 0))
        }
        return differences

    def check_specific_target_loss(self, diff_comp: dict, target_losses: list) -> bool:
        for target_loss in target_losses:
            target_elements = dict(mass.Composition(target_loss))
            print(f"Checking target loss - Differences: {diff_comp}, Target: {target_elements}")
            if all(diff_comp.get(k, 0) == target_elements.get(k, 0) for k in target_elements):
                return True
        return False

    def clean_title(self, title):
        return str(title).strip()

    def process_files(self, input_path_csv, input_path_mgf, output_path_csv, db_path,
                    target_loss1, target_loss2, target_loss3, target_loss4, target_loss5, target_loss6,
                    charge, ppm):
        df = pd.read_csv(input_path_csv)
        df.set_index('row ID', inplace=True)
        df.index = df.index.map(self.clean_title)
        formula_ranges = self.read_formulas_from_NLdb(db_path, charge, ppm)

        lossname1 = self.entry_lossname1.get()
        lossname2 = self.entry_lossname2.get()
        lossname3 = self.entry_lossname3.get()
        lossname4 = self.entry_lossname4.get()
        lossname5 = self.entry_lossname5.get()
        lossname6 = self.entry_lossname6.get()

        df[f'{lossname1}'] = 'none'
        df[f'{lossname2}'] = 'none'
        df[f'{lossname3}'] = 'none'
        df[f'{lossname4}'] = 'none'
        df[f'{lossname5}'] = 'none'
        df[f'{lossname6}'] = 'none'

        spectra_dict = {}
        with mgf.read(input_path_mgf) as spectra:
            for spectrum in spectra:
                cleaned_title = self.clean_title(spectrum['params']['title'])
                spectra_dict[cleaned_title] = spectrum

        total_rows = len(df.index)
        for i, row_id in enumerate(df.index):
            if row_id in spectra_dict:
                spectrum = spectra_dict[row_id]
                precursor_ion = spectrum['params']['pepmass'][0]
                df.at[row_id, 'row m/z'] = precursor_ion
                product_ions = spectrum['m/z array']
                formulas = [self.find_matching_NLformulas(mz, formula_ranges) for mz in product_ions]

                # loss 1
                carbon_range1_1_max = self.entry_carbonRange1_1_max.get()
                carbon_range1_1_min = self.entry_carbonRange1_1_min.get()
                carbon_range1_2_max = self.entry_carbonRange1_2_max.get()
                carbon_range1_2_min = self.entry_carbonRange1_2_min.get()

                # loss 2
                carbon_range2_1_max = self.entry_carbonRange2_1_max.get()
                carbon_range2_1_min = self.entry_carbonRange2_1_min.get()
                carbon_range2_2_max = self.entry_carbonRange2_2_max.get()
                carbon_range2_2_min = self.entry_carbonRange2_2_min.get()

                # loss 3
                carbon_range3_1_max = self.entry_carbonRange3_1_max.get()
                carbon_range3_1_min = self.entry_carbonRange3_1_min.get()
                carbon_range3_2_max = self.entry_carbonRange3_2_max.get()
                carbon_range3_2_min = self.entry_carbonRange3_2_min.get()

                # loss 4
                carbon_range4_1_max = self.entry_carbonRange4_1_max.get()
                carbon_range4_1_min = self.entry_carbonRange4_1_min.get()
                carbon_range4_2_max = self.entry_carbonRange4_2_max.get()
                carbon_range4_2_min = self.entry_carbonRange4_2_min.get()

                # loss 5
                carbon_range5_1_max = self.entry_carbonRange5_1_max.get()
                carbon_range5_1_min = self.entry_carbonRange5_1_min.get()
                carbon_range5_2_max = self.entry_carbonRange5_2_max.get()
                carbon_range5_2_min = self.entry_carbonRange5_2_min.get()

                # loss 6
                carbon_range6_1_max = self.entry_carbonRange6_1_max.get()
                carbon_range6_1_min = self.entry_carbonRange6_1_min.get()
                carbon_range6_2_max = self.entry_carbonRange6_2_max.get()
                carbon_range6_2_min = self.entry_carbonRange6_2_min.get()

                self.progress_bar.set((i + 1) / total_rows)

                # loss 1
                carbon_range1_1_max_formulas = [f for f in formulas if
                                                f != 'none' and f.startswith(f'{carbon_range1_1_max}')]
                carbon_range1_1_min_formulas = [f for f in formulas if
                                                f != 'none' and f.startswith(f'{carbon_range1_1_min}')]
                carbon_range1_2_max_formulas = [f for f in formulas if
                                                f != 'none' and f.startswith(f'{carbon_range1_2_max}')]
                carbon_range1_2_min_formulas = [f for f in formulas if
                                                f != 'none' and f.startswith(f'{carbon_range1_2_min}')]

                # loss 2
                carbon_range2_1_max_formulas = [f for f in formulas if
                                                f != 'none' and f.startswith(f'{carbon_range2_1_max}')]
                carbon_range2_1_min_formulas = [f for f in formulas if
                                                f != 'none' and f.startswith(f'{carbon_range2_1_min}')]
                carbon_range2_2_max_formulas = [f for f in formulas if
                                                f != 'none' and f.startswith(f'{carbon_range2_2_max}')]
                carbon_range2_2_min_formulas = [f for f in formulas if
                                                f != 'none' and f.startswith(f'{carbon_range2_2_min}')]

                # loss 3
                carbon_range3_1_max_formulas = [f for f in formulas if
                                                f != 'none' and f.startswith(f'{carbon_range3_1_max}')]
                carbon_range3_1_min_formulas = [f for f in formulas if
                                                f != 'none' and f.startswith(f'{carbon_range3_1_min}')]
                carbon_range3_2_max_formulas = [f for f in formulas if
                                                f != 'none' and f.startswith(f'{carbon_range3_2_max}')]
                carbon_range3_2_min_formulas = [f for f in formulas if
                                                f != 'none' and f.startswith(f'{carbon_range3_2_min}')]

                # loss 4
                carbon_range4_1_max_formulas = [f for f in formulas if
                                                f != 'none' and f.startswith(f'{carbon_range4_1_max}')]
                carbon_range4_1_min_formulas = [f for f in formulas if
                                                f != 'none' and f.startswith(f'{carbon_range4_1_min}')]
                carbon_range4_2_max_formulas = [f for f in formulas if
                                                f != 'none' and f.startswith(f'{carbon_range4_2_max}')]
                carbon_range4_2_min_formulas = [f for f in formulas if
                                                f != 'none' and f.startswith(f'{carbon_range4_2_min}')]

                # loss 5
                carbon_range5_1_max_formulas = [f for f in formulas if
                                                f != 'none' and f.startswith(f'{carbon_range5_1_max}')]
                carbon_range5_1_min_formulas = [f for f in formulas if
                                                f != 'none' and f.startswith(f'{carbon_range5_1_min}')]
                carbon_range5_2_max_formulas = [f for f in formulas if
                                                f != 'none' and f.startswith(f'{carbon_range5_2_max}')]
                carbon_range5_2_min_formulas = [f for f in formulas if
                                                f != 'none' and f.startswith(f'{carbon_range5_2_min}')]

                # loss 6
                carbon_range6_1_max_formulas = [f for f in formulas if
                                                f != 'none' and f.startswith(f'{carbon_range6_1_max}')]
                carbon_range6_1_min_formulas = [f for f in formulas if
                                                f != 'none' and f.startswith(f'{carbon_range6_1_min}')]
                carbon_range6_2_max_formulas = [f for f in formulas if
                                                f != 'none' and f.startswith(f'{carbon_range6_2_max}')]
                carbon_range6_2_min_formulas = [f for f in formulas if
                                                f != 'none' and f.startswith(f'{carbon_range6_2_min}')]

                # Print the identified formulas for debugging and verification.
                print(f"Spectrum Title: {row_id}, Precursor Ion m/z: {precursor_ion}")
                # loss 1
                print(f"Identified {carbon_range1_1_max} Formulas: {carbon_range1_1_max_formulas}")
                print(f"Identified {carbon_range1_1_min} Formulas: {carbon_range1_1_min_formulas}")
                print(f"Identified {carbon_range1_2_max} Formulas: {carbon_range1_2_max_formulas}")
                print(f"Identified {carbon_range1_2_min} Formulas: {carbon_range1_2_min_formulas}")
                # loss 2
                print(f"Identified {carbon_range2_1_max} Formulas: {carbon_range2_1_max_formulas}")
                print(f"Identified {carbon_range2_1_min} Formulas: {carbon_range2_1_min_formulas}")
                print(f"Identified {carbon_range2_2_max} Formulas: {carbon_range2_2_max_formulas}")
                print(f"Identified {carbon_range2_2_min} Formulas: {carbon_range2_2_min_formulas}")
                # loss 3
                print(f"Identified {carbon_range3_1_max} Formulas: {carbon_range3_1_max_formulas}")
                print(f"Identified {carbon_range3_1_min} Formulas: {carbon_range3_1_min_formulas}")
                print(f"Identified {carbon_range3_2_max} Formulas: {carbon_range3_2_max_formulas}")
                print(f"Identified {carbon_range3_2_min} Formulas: {carbon_range3_2_min_formulas}")
                # loss 4
                print(f"Identified {carbon_range4_1_max} Formulas: {carbon_range4_1_max_formulas}")
                print(f"Identified {carbon_range4_1_min} Formulas: {carbon_range4_1_min_formulas}")
                print(f"Identified {carbon_range4_2_max} Formulas: {carbon_range4_2_max_formulas}")
                print(f"Identified {carbon_range4_2_min} Formulas: {carbon_range4_2_min_formulas}")
                # loss 5
                print(f"Identified {carbon_range5_1_max} Formulas: {carbon_range5_1_max_formulas}")
                print(f"Identified {carbon_range5_1_min} Formulas: {carbon_range5_1_min_formulas}")
                print(f"Identified {carbon_range5_2_max} Formulas: {carbon_range5_2_max_formulas}")
                print(f"Identified {carbon_range5_2_min} Formulas: {carbon_range5_2_min_formulas}")
                # loss 6
                print(f"Identified {carbon_range6_1_max} Formulas: {carbon_range6_1_max_formulas}")
                print(f"Identified {carbon_range6_1_min} Formulas: {carbon_range6_1_min_formulas}")
                print(f"Identified {carbon_range6_2_max} Formulas: {carbon_range6_2_max_formulas}")
                print(f"Identified {carbon_range6_2_min} Formulas: {carbon_range6_2_min_formulas}")

                # loss 1
                for loss1 in target_loss1:
                    for carbon_range1_1_max_formula in carbon_range1_1_max_formulas:
                        for carbon_range1_1_min_formula in carbon_range1_1_min_formulas:
                            diff = self.specific_formula_difference(carbon_range1_1_max_formula,
                                                                    carbon_range1_1_min_formula)
                            print(
                                f"Checking {carbon_range1_1_max} formula {carbon_range1_1_max_formula} against {carbon_range1_1_min} formula {carbon_range1_1_min_formula}: Difference {diff}")
                            if self.check_specific_target_loss(diff, [loss1]):
                                print(
                                    f"Product ions pair with target loss {loss1}: {carbon_range1_1_max_formula} and {carbon_range1_1_min_formula}")
                                df.at[row_id, f'{lossname1}'] = loss1

                for loss1 in target_loss1:
                    for carbon_range1_2_max_formula in carbon_range1_2_max_formulas:
                        for carbon_range1_2_min_formula in carbon_range1_2_min_formulas:
                            diff = self.specific_formula_difference(carbon_range1_2_max_formula,
                                                                    carbon_range1_2_min_formula)
                            print(
                                f"Checking {carbon_range1_2_max} formula {carbon_range1_2_max_formula} against {carbon_range1_2_min} formula {carbon_range1_2_min_formula}: Difference {diff}")
                            if self.check_specific_target_loss(diff, [loss1]):
                                print(
                                    f"Product ions pair with target loss {loss1}: {carbon_range1_2_max_formula} and {carbon_range1_2_min_formula}")
                                df.at[row_id, f'{lossname1}'] = loss1

                # loss 2
                for loss2 in target_loss2:
                    for carbon_range2_1_max_formula in carbon_range2_1_max_formulas:
                        for carbon_range2_1_min_formula in carbon_range2_1_min_formulas:
                            diff = self.specific_formula_difference(carbon_range2_1_max_formula,
                                                                    carbon_range2_1_min_formula)
                            print(
                                f"Checking {carbon_range2_1_max} formula {carbon_range2_1_max_formula} against {carbon_range2_1_min} formula {carbon_range2_1_min_formula}: Difference {diff}")
                            if self.check_specific_target_loss(diff, [loss2]):
                                print(
                                    f"Product ions pair with target loss {loss2}: {carbon_range2_1_max_formula} and {carbon_range2_1_min_formula}")
                                df.at[row_id, f'{lossname2}'] = loss2

                for loss2 in target_loss2:
                    for carbon_range2_2_max_formula in carbon_range2_2_max_formulas:
                        for carbon_range2_2_min_formula in carbon_range2_2_min_formulas:
                            diff = self.specific_formula_difference(carbon_range2_2_max_formula,
                                                                    carbon_range2_2_min_formula)
                            print(
                                f"Checking {carbon_range2_2_max} formula {carbon_range2_2_max_formula} against {carbon_range2_2_min} formula {carbon_range2_2_min_formula}: Difference {diff}")
                            if self.check_specific_target_loss(diff, [loss2]):
                                print(
                                    f"Product ions pair with target loss {loss2}: {carbon_range2_2_max_formula} and {carbon_range2_2_min_formula}")
                                df.at[row_id, f'{lossname2}'] = loss2

                # loss 3
                for loss3 in target_loss3:
                    for carbon_range3_1_max_formula in carbon_range3_1_max_formulas:
                        for carbon_range3_1_min_formula in carbon_range3_1_min_formulas:
                            diff = self.specific_formula_difference(carbon_range3_1_max_formula,
                                                                    carbon_range3_1_min_formula)
                            print(
                                f"Checking {carbon_range3_1_max} formula {carbon_range3_1_max_formula} against {carbon_range3_1_min} formula {carbon_range3_1_min_formula}: Difference {diff}")
                            if self.check_specific_target_loss(diff, [loss3]):
                                print(
                                    f"Product ions pair with target loss {loss3}: {carbon_range3_1_max_formula} and {carbon_range3_1_min_formula}")
                                df.at[row_id, f'{lossname3}'] = loss3

                for loss3 in target_loss3:
                    for carbon_range3_2_max_formula in carbon_range3_2_max_formulas:
                        for carbon_range3_2_min_formula in carbon_range3_2_min_formulas:
                            diff = self.specific_formula_difference(carbon_range3_2_max_formula,
                                                                    carbon_range3_2_min_formula)
                            print(
                                f"Checking {carbon_range3_2_max} formula {carbon_range3_2_max_formula} against {carbon_range3_2_min} formula {carbon_range3_2_min_formula}: Difference {diff}")
                            if self.check_specific_target_loss(diff, [loss3]):
                                print(
                                    f"Product ions pair with target loss {loss3}: {carbon_range3_2_max_formula} and {carbon_range3_2_min_formula}")
                                df.at[row_id, f'{lossname3}'] = loss3

                # loss 4
                for loss4 in target_loss4:
                    for carbon_range4_1_max_formula in carbon_range4_1_max_formulas:
                        for carbon_range4_1_min_formula in carbon_range4_1_min_formulas:
                            diff = self.specific_formula_difference(carbon_range4_1_max_formula,
                                                                    carbon_range4_1_min_formula)
                            print(
                                f"Checking {carbon_range4_1_max} formula {carbon_range4_1_max_formula} against {carbon_range4_1_min} formula {carbon_range4_1_min_formula}: Difference {diff}")
                            if self.check_specific_target_loss(diff, [loss4]):
                                print(
                                    f"Product ions pair with target loss {loss4}: {carbon_range4_1_max_formula} and {carbon_range4_1_min_formula}")
                                df.at[row_id, f'{lossname4}'] = loss4

                for loss4 in target_loss4:
                    for carbon_range4_2_max_formula in carbon_range4_2_max_formulas:
                        for carbon_range4_2_min_formula in carbon_range4_2_min_formulas:
                            diff = self.specific_formula_difference(carbon_range4_2_max_formula,
                                                                    carbon_range4_2_min_formula)
                            print(
                                f"Checking {carbon_range4_2_max} formula {carbon_range4_2_max_formula} against {carbon_range4_2_min} formula {carbon_range4_2_min_formula}: Difference {diff}")
                            if self.check_specific_target_loss(diff, [loss4]):
                                print(
                                    f"Product ions pair with target loss {loss4}: {carbon_range4_2_max_formula} and {carbon_range4_2_min_formula}")
                                df.at[row_id, f'{lossname4}'] = loss4

                # loss 5
                for loss5 in target_loss5:
                    for carbon_range5_1_max_formula in carbon_range5_1_max_formulas:
                        for carbon_range5_1_min_formula in carbon_range5_1_min_formulas:
                            diff = self.specific_formula_difference(carbon_range5_1_max_formula,
                                                                    carbon_range5_1_min_formula)
                            print(
                                f"Checking {carbon_range5_1_max} formula {carbon_range5_1_max_formula} against {carbon_range5_1_min} formula {carbon_range5_1_min_formula}: Difference {diff}")
                            if self.check_specific_target_loss(diff, [loss5]):
                                print(
                                    f"Product ions pair with target loss {loss5}: {carbon_range5_1_max_formula} and {carbon_range5_1_min_formula}")
                                df.at[row_id, f'{lossname5}'] = loss5

                for loss5 in target_loss5:
                    for carbon_range5_2_max_formula in carbon_range5_2_max_formulas:
                        for carbon_range5_2_min_formula in carbon_range5_2_min_formulas:
                            diff = self.specific_formula_difference(carbon_range5_2_max_formula,
                                                                    carbon_range5_2_min_formula)
                            print(
                                f"Checking {carbon_range5_2_max} formula {carbon_range5_2_max_formula} against {carbon_range5_2_min} formula {carbon_range5_2_min_formula}: Difference {diff}")
                            if self.check_specific_target_loss(diff, [loss5]):
                                print(
                                    f"Product ions pair with target loss {loss5}: {carbon_range5_2_max_formula} and {carbon_range5_2_min_formula}")
                                df.at[row_id, f'{lossname5}'] = loss5

                # loss 6
                for loss6 in target_loss6:
                    for carbon_range6_1_max_formula in carbon_range6_1_max_formulas:
                        for carbon_range6_1_min_formula in carbon_range6_1_min_formulas:
                            diff = self.specific_formula_difference(carbon_range6_1_max_formula,
                                                                    carbon_range6_1_min_formula)
                            print(
                                f"Checking {carbon_range6_1_max} formula {carbon_range6_1_max_formula} against {carbon_range6_1_min} formula {carbon_range6_1_min_formula}: Difference {diff}")
                            if self.check_specific_target_loss(diff, [loss6]):
                                print(
                                    f"Product ions pair with target loss {loss6}: {carbon_range6_1_max_formula} and {carbon_range6_1_min_formula}")
                                df.at[row_id, f'{lossname6}'] = loss6

                for loss6 in target_loss6:
                    for carbon_range6_2_max_formula in carbon_range6_2_max_formulas:
                        for carbon_range6_2_min_formula in carbon_range6_2_min_formulas:
                            diff = self.specific_formula_difference(carbon_range6_2_max_formula,
                                                                    carbon_range6_2_min_formula)
                            print(
                                f"Checking {carbon_range6_2_max} formula {carbon_range6_2_max_formula} against {carbon_range6_2_min} formula {carbon_range6_2_min_formula}: Difference {diff}")
                            if self.check_specific_target_loss(diff, [loss6]):
                                print(
                                    f"Product ions pair with target loss {loss6}: {carbon_range6_2_max_formula} and {carbon_range6_2_min_formula}")
                                df.at[row_id, f'{lossname6}'] = loss6

        def composition_to_formula(composition_str):
            clean_str = composition_str.replace("Composition(", "").replace(")", "").replace("'", "")
            matches = re.findall(r'([A-Z][a-z]*): (\d+)', clean_str)
            formula = ''.join([f"{elem}{count}" for elem, count in matches])
            return formula

        comp_index = df.columns.get_loc('same composition') + 1
        df.insert(comp_index, f'{lossname1}', df.pop(f'{lossname1}'))
        df.insert(comp_index + 1, f'{lossname2}', df.pop(f'{lossname2}'))
        df.insert(comp_index + 2, f'{lossname3}', df.pop(f'{lossname3}'))
        df.insert(comp_index + 3, f'{lossname4}', df.pop(f'{lossname4}'))
        df.insert(comp_index + 4, f'{lossname5}', df.pop(f'{lossname5}'))
        df.insert(comp_index + 4, f'{lossname6}', df.pop(f'{lossname6}'))

        df[f'{lossname1}'] = df[f'{lossname1}'].apply(lambda x: composition_to_formula(str(x)) if x != 'none' else x)
        df[f'{lossname2}'] = df[f'{lossname2}'].apply(lambda x: composition_to_formula(str(x)) if x != 'none' else x)
        df[f'{lossname3}'] = df[f'{lossname3}'].apply(lambda x: composition_to_formula(str(x)) if x != 'none' else x)
        df[f'{lossname4}'] = df[f'{lossname4}'].apply(lambda x: composition_to_formula(str(x)) if x != 'none' else x)
        df[f'{lossname5}'] = df[f'{lossname5}'].apply(lambda x: composition_to_formula(str(x)) if x != 'none' else x)
        df[f'{lossname6}'] = df[f'{lossname6}'].apply(lambda x: composition_to_formula(str(x)) if x != 'none' else x)

        df.reset_index().to_csv(output_path_csv, index=False)
        self.progress_bar.set(1.0)


    '''feature_ion_extractor_tab'''
    def create_feature_ion_extractor_tab(self):
        feature_ion_extractor_tab = self.tab("5.Feature Ion Extract")

        # Create a canvas
        canvas = tk.Canvas(feature_ion_extractor_tab)
        canvas.pack(side="left", fill="both", expand=True)

        # Add a scrollbar
        scrollbar = tk.Scrollbar(feature_ion_extractor_tab, orient="vertical", command=canvas.yview)
        scrollbar.pack(side="right", fill="y")

        # Configure canvas with scrollbar
        canvas.configure(yscrollcommand=scrollbar.set)
        canvas.bind("<Configure>", lambda e: canvas.configure(scrollregion=canvas.bbox("all")))

        # Create a frame inside the canvas
        self.frame_inside_canvas = ctk.CTkFrame(canvas)
        canvas.create_window((0, 0), window=self.frame_inside_canvas, anchor='nw')

        # CTkFrame of Import MGF
        frame_ImportMGF = ctk.CTkFrame(self.frame_inside_canvas, fg_color="transparent")
        frame_ImportMGF.pack(fill='x', padx=10, pady=5)
        # Section for MGF file input
        label_extractor_mgf = ctk.CTkLabel(frame_ImportMGF, text="Import query MS2 Spectra:", font=('Arial', 14, 'bold'))
        label_extractor_mgf.pack(side='left', padx=(0, 10))
        self.entry_FIextractor_mgf_path = ctk.CTkEntry(frame_ImportMGF, placeholder_text="Input MGF file path", width=490)
        self.entry_FIextractor_mgf_path.pack(side='left', padx=(0, 5))
        self.button_FIextractor_mgf_browse = ctk.CTkButton(frame_ImportMGF, text="Browse file", command=self.browse_FIextractor_mgf, fg_color="#6598D3", hover_color="#5578B3")
        self.button_FIextractor_mgf_browse.pack(padx=10, pady=5, anchor='e')

        # CTkFrame of Import CSV
        frame_ImportCSV = ctk.CTkFrame(self.frame_inside_canvas, fg_color="transparent")
        frame_ImportCSV.pack(fill='x', padx=10, pady=5)
        # Section for CSV file input
        label_extractor_csv = ctk.CTkLabel(frame_ImportCSV, text="Import query MS Spectra:", font=('Arial', 14, 'bold'))
        label_extractor_csv.pack(side='left', padx=(0, 10))
        self.entry_FIextractor_csv_path = ctk.CTkEntry(frame_ImportCSV, placeholder_text="Input CSV file path", width=500)
        self.entry_FIextractor_csv_path.pack(side='left', padx=(0, 5))
        self.button_FIextractor_csv_browse = ctk.CTkButton(frame_ImportCSV, text="Browse file", command=self.browse_FIextractor_csv, fg_color="#6598D3", hover_color="#5578B3")
        self.button_FIextractor_csv_browse.pack(padx=10, pady=5, anchor='e')

        # CTkFrame of Output
        frame_Output = ctk.CTkFrame(self.frame_inside_canvas, fg_color="transparent")
        frame_Output.pack(fill='x', padx=10, pady=5)
        label_output = ctk.CTkLabel(frame_Output, text="Export Result Files:", font=('Arial', 14, 'bold'))
        label_output.pack(side='left', padx=(0, 10))
        self.entry_FIextractor_output_path = ctk.CTkEntry(frame_Output, placeholder_text="Output csv file path, e.g. C:/Users/xxx/", width=330)
        self.entry_FIextractor_output_path.pack(side='left', padx=(0, 5))
        self.entry_FIextractor_output_csv_filename = ctk.CTkEntry(frame_Output, placeholder_text="Output csv file name, e.g. xxx.csv", width=205)
        self.entry_FIextractor_output_csv_filename.pack(side='left', padx=(0, 5))
        self.button_FIextractor_calculate = ctk.CTkButton(frame_Output, text="Extract Feature Ion", command=self.FIextract, fg_color="#6598D3", hover_color="#5578B3")
        self.button_FIextractor_calculate.pack(padx=10, pady=5, anchor='e')

        # Create initial ion set inside the canvas
        self.add_ion_set()

        # Adding + button to add new ion sets
        add_button = ctk.CTkButton(self.frame_inside_canvas, text="+ Add Specific Ion Set", command=self.add_ion_set, fg_color="#6598D3", hover_color="#5578B3")
        add_button.pack(side='bottom', pady=10)

    def add_ion_set(self):
        ion_set_frame = ctk.CTkFrame(self.frame_inside_canvas)
        ion_set_frame.pack(fill='x', padx=10, pady=5)

        ion_set_frame.entries = {}

        set_number = len(self.specific_ion_sets) + 1
        self.specific_ion_sets.append(ion_set_frame)

        Specific_Ion_Set_label = ctk.CTkLabel(ion_set_frame, text=f"Specific Ion Set {set_number}:",
                                            font=('Arial', 14, 'bold'))
        Specific_Ion_Set_label.pack(side='left', padx=(0, 10))

        # Add entry for formulas
        ion_set_entry = ctk.CTkEntry(ion_set_frame,
                                        placeholder_text=f"Specific Ion Set {set_number}, e.g.,C10H16O1, C10H14O1",
                                        width=460)
        ion_set_entry.pack(side='left', padx=(0, 5))

        # Add entry for name
        ion_set_name_entry = ctk.CTkEntry(ion_set_frame, placeholder_text=f"Name {set_number}, e.g.,C10 Ion", width=150)
        ion_set_name_entry.pack(side='left', padx=(0, 5))

        # Store formulas and name in entries
        ion_set_frame.entries['formulas'] = ion_set_entry
        ion_set_frame.entries['name'] = ion_set_name_entry

        # Add parameters frame
        self.add_ion_set_parameters(ion_set_frame)

    def add_ion_set_parameters(self, ion_set_frame):
        frame_Specific_Ion_Set_para1 = ctk.CTkFrame(ion_set_frame, fg_color="transparent")
        frame_Specific_Ion_Set_para1.pack(fill='x', padx=10, pady=5)

        label_massrange_min = ctk.CTkLabel(frame_Specific_Ion_Set_para1, text="・Mass Range of Product Ions (m/z):",
                                            font=('Arial', 14))
        label_massrange_min.pack(side='left', padx=(0, 10))
        entry_mass_range_min_Set = ctk.CTkEntry(frame_Specific_Ion_Set_para1, placeholder_text="min, e.g.40", width=203)
        entry_mass_range_min_Set.pack(side='left', padx=(0, 5))

        label_massrange_max = ctk.CTkLabel(frame_Specific_Ion_Set_para1, text="to", font=('Arial', 14))
        label_massrange_max.pack(side='left', padx=(0, 5))

        max_mass_option = ctk.CTkOptionMenu(frame_Specific_Ion_Set_para1, values=["Use Precursor Ion", "Specify Value"],
                                            fg_color="#F5F5F5", button_color="#BEBEBE", text_color="#333333", width=160)
        max_mass_option.pack(side='left', padx=(0, 5))
        entry_mass_range_max_Set = ctk.CTkEntry(frame_Specific_Ion_Set_para1, placeholder_text="max, e.g.250",
                                                width=190)
        entry_mass_range_max_Set.pack(side='left', padx=(0, 5))

        frame_Specific_Ion_Set_para2 = ctk.CTkFrame(ion_set_frame, fg_color="transparent")
        frame_Specific_Ion_Set_para2.pack(fill='x', padx=10, pady=5)

        label_charge_Set = ctk.CTkLabel(frame_Specific_Ion_Set_para2, text="・Charge:", font=('Arial', 14))
        label_charge_Set.pack(side='left', padx=(0, 5))
        entry_charge_Set = ctk.CTkEntry(frame_Specific_Ion_Set_para2, placeholder_text="eg.+1/-1", width=76)
        entry_charge_Set.pack(side='left', padx=(0, 5))

        label_massTolerance_Set = ctk.CTkLabel(frame_Specific_Ion_Set_para2, text="・Mass Tolerence (ppm):",
                                            font=('Arial', 14))
        label_massTolerance_Set.pack(side='left', padx=(0, 5))
        entry_mass_tolerance_Set = ctk.CTkEntry(frame_Specific_Ion_Set_para2, placeholder_text="e.g., 5", width=76)
        entry_mass_tolerance_Set.pack(side='left', padx=(0, 5))

        label_IntensityNor_Set = ctk.CTkLabel(frame_Specific_Ion_Set_para2, text="・Intensity Normalization:",
                                            font=('Arial', 14))
        label_IntensityNor_Set.pack(side='left', padx=(0, 5))
        entry_IntensityNor_Set = ctk.CTkEntry(frame_Specific_Ion_Set_para2, placeholder_text="e.g., 0.1", width=76)
        entry_IntensityNor_Set.pack(side='left', padx=(0, 5))

        label_TopNumber_Set = ctk.CTkLabel(frame_Specific_Ion_Set_para2, text="・Top Number:", font=('Arial', 14))
        label_TopNumber_Set.pack(side='left', padx=(0, 5))
        entry_TopNumber_type = ctk.CTkEntry(frame_Specific_Ion_Set_para2, placeholder_text="e.g., 1, 2, 3", width=76)
        entry_TopNumber_type.pack(side='left', padx=(0, 5))

        # Update the entries dictionary with these new entries
        ion_set_frame.entries.update({
            'mass_min': entry_mass_range_min_Set,
            'mass_max': entry_mass_range_max_Set,
            'charge': entry_charge_Set,
            'tolerance': entry_mass_tolerance_Set,
            'intensity_normalization': entry_IntensityNor_Set,
            'top_number': entry_TopNumber_type,
            'max_mass_option': max_mass_option
        })

    def browse_FIextractor_mgf(self):
        path = filedialog.askopenfilename(filetypes=[("MGF files", "*.mgf")])
        if path:
            self.entry_FIextractor_mgf_path.delete(0, tk.END)
            self.entry_FIextractor_mgf_path.insert(0, path)

    def browse_FIextractor_csv(self):
        path = filedialog.askopenfilename(filetypes=[("CSV files", "*.csv")])
        if path:
            self.entry_FIextractor_csv_path.delete(0, tk.END)
            self.entry_FIextractor_csv_path.insert(0, path)

    def FIextract(self):
        try:
            output_dir = self.entry_FIextractor_output_path.get().strip()
            output_filename = self.entry_FIextractor_output_csv_filename.get().strip()
            full_output_path = os.path.join(output_dir, output_filename)

            if not os.path.exists(output_dir):
                os.makedirs(output_dir)

            threading.Thread(target=self.process_FeatureIonExtract_parallel, args=(
                self.entry_FIextractor_csv_path.get(),
                self.entry_FIextractor_mgf_path.get(),
                full_output_path
            )).start()
            messagebox.showinfo("Process Complete", "Feature ions have been successfully extracted.")
        except Exception as e:
            messagebox.showerror("Error", f"Failed to save file: {str(e)}")

    def detect_formulas_in_spectra(self, spectrum, specific_formulas, params):
        matched_formulas = {}
        try:
            min_mass = int(params['mass_min'].get())
            charge = int(params['charge'].get())
            tolerance_ppm = int(params['tolerance'].get())
            norScore = float(params['intensity_normalization'].get())
            top_number = int(params['top_number'].get())
            max_mass_option = params['max_mass_option'].get()

            title = spectrum['params']['title']
            if max_mass_option == "Specify Value":
                max_mass = int(params['mass_max'].get())
            else:
                max_mass = spectrum['params']['pepmass'][0]

            product_ions = spectrum['m/z array']
            intensity = spectrum['intensity array']
            ind = np.nonzero((product_ions > min_mass) & (product_ions < max_mass) & (intensity > norScore))[0]

            if len(ind) > 0:
                intensity_max = np.max(intensity[ind])
                intensity_min = np.min(intensity[ind])
                normalized_intensity = np.zeros_like(intensity)
                normalized_intensity[ind] = (intensity[ind] - intensity_min) / (intensity_max - intensity_min)
                formula_intensity_pairs = []

                for formula in specific_formulas:
                    mz = mass.calculate_mass(formula=formula, charge=charge)
                    mz_min = mz - mz * tolerance_ppm / 1e6
                    mz_max = mz + mz * tolerance_ppm / 1e6

                    for i in ind:
                        if mz_min <= product_ions[i] <= mz_max and normalized_intensity[i] > norScore:
                            formula_intensity_pairs.append((formula, normalized_intensity[i]))


                formula_intensity_pairs.sort(key=lambda x: x[1], reverse=True)
                top_formulas = [pair[0] for pair in formula_intensity_pairs[:top_number]]
                if top_formulas:
                    matched_formulas[int(title)] = ', '.join(top_formulas)

        except Exception as e:
            print("Error during formula detection:", str(e))

        return matched_formulas

    def process_FeatureIonExtract_parallel(self, path_csv, path_mgf, full_output_path):
        try:
            df = pd.read_csv(path_csv)
            df.set_index('row ID', inplace=True)

            spectra_dict = {}
            with mgf.read(path_mgf) as spectra:
                for spectrum in spectra:
                    title = int(spectrum['params']['title'])
                    spectra_dict[title] = spectrum

            specific_formulas = {}
            for ion_set_frame in self.specific_ion_sets:
                ion_set_name = ion_set_frame.entries['name'].get()
                specific_formulas[ion_set_name] = [x.strip() for x in ion_set_frame.entries['formulas'].get().split(',')]

            def process_ion_set(ion_set_frame, ion_set_name):
                ion_set_params = ion_set_frame.entries
                matched_formulas = {}
                for title, spectrum in spectra_dict.items():
                    result = self.detect_formulas_in_spectra(spectrum, specific_formulas[ion_set_name], ion_set_params)
                    if result:
                        matched_formulas.update(result)
                return ion_set_name, matched_formulas

            with ThreadPoolExecutor() as executor:
                futures = [executor.submit(process_ion_set, ion_set_frame, ion_set_name) for ion_set_frame, ion_set_name in zip(self.specific_ion_sets, specific_formulas.keys())]

                for future in futures:
                    ion_set_name, matched_formulas = future.result()
                    df[ion_set_name] = df.index.map(lambda x: matched_formulas.get(x, 'none'))

            df.reset_index(inplace=True)
            df.to_csv(full_output_path, index=False)
            print("Data processed and saved to", full_output_path)
        except Exception as e:
            print("Error during processing:", str(e))
            raise


    '''cosine_score_calculator_tab'''
    def create_cosine_score_calculator_tab(self):
        cosine_score_calculator_tab = self.tab("6.Cosine Score Calculate")

        # Create a canvas
        canvas = tk.Canvas(cosine_score_calculator_tab)
        canvas.pack(side="left", fill="both", expand=True)

        # Create a frame inside the canvas
        frame_inside_canvas = ctk.CTkFrame(canvas)
        canvas.create_window((0, 0), window=frame_inside_canvas, anchor='nw')

        # CTkFrame of Import MGF
        frame_ImportCosine_query_mgf = ctk.CTkFrame(frame_inside_canvas, fg_color="transparent")
        frame_ImportCosine_query_mgf.pack(fill='x', padx=10, pady=5)
        # Section for Query MGF file input
        label_cosine_query_mgf = ctk.CTkLabel(frame_ImportCosine_query_mgf, text="Import Query MS2 Spectra (MGF file):", font=('Arial', 14, 'bold'))
        label_cosine_query_mgf.pack(side='left', padx=(0, 10))
        self.entry_cosine_query_mgf_path = ctk.CTkEntry(frame_ImportCosine_query_mgf, placeholder_text="Input query MGF file path", width=400)
        self.entry_cosine_query_mgf_path.pack(side='left', padx=(0, 5))
        self.button_cosine_query_mgf_browse = ctk.CTkButton(frame_ImportCosine_query_mgf, text="Browse file",
                                                    command=self.browse_query_mgf, fg_color="#6598D3", hover_color="#5578B3")
        self.button_cosine_query_mgf_browse.pack(padx=10, pady=5, anchor='e')

        # CTkFrame of Import CSV
        frame_ImportCosine_query_csv = ctk.CTkFrame(frame_inside_canvas, fg_color="transparent")
        frame_ImportCosine_query_csv.pack(fill='x', padx=10, pady=5)
        # Section for Query CSV file input
        label_cosine_query_csv = ctk.CTkLabel(frame_ImportCosine_query_csv, text="Import Query MS Spectra (CSV file):", font=('Arial', 14, 'bold'))
        label_cosine_query_csv.pack(side='left', padx=(0, 10))
        self.entry_cosine_query_csv_path = ctk.CTkEntry(frame_ImportCosine_query_csv, placeholder_text="Input query CSV file path", width=410)
        self.entry_cosine_query_csv_path.pack(side='left', padx=(0, 5))
        self.button_cosine_query_csv_browse = ctk.CTkButton(frame_ImportCosine_query_csv, text="Browse file",
                                            command=self.browse_query_csv, fg_color="#6598D3", hover_color="#5578B3")
        self.button_cosine_query_csv_browse.pack(padx=10, pady=5, anchor='e')

        # CTkFrame of Import MGF
        frame_ImportCosine_std_mgf = ctk.CTkFrame(frame_inside_canvas, fg_color="transparent")
        frame_ImportCosine_std_mgf.pack(fill='x', padx=10, pady=5)
        # Section for Standards MGF file input
        label_cosine_std_mgf = ctk.CTkLabel(frame_ImportCosine_std_mgf, text="Import Standards MS2 Spectra (MGF file):", font=('Arial', 14, 'bold'))
        label_cosine_std_mgf.pack(side='left', padx=(0, 10))
        self.entry_cosine_std_mgf_path = ctk.CTkEntry(frame_ImportCosine_std_mgf, placeholder_text="Input Standards MGF file path", width=375)
        self.entry_cosine_std_mgf_path.pack(side='left', padx=(0, 5))
        self.button_cosine_std_mgf_browse = ctk.CTkButton(frame_ImportCosine_std_mgf, text="Browse file",
                                                    command=self.browse_std_mgf, fg_color="#6598D3", hover_color="#5578B3")
        self.button_cosine_std_mgf_browse.pack(padx=10, pady=5, anchor='e')

        # CTkFrame of Import CSV
        frame_ImportCosine_std_csv = ctk.CTkFrame(frame_inside_canvas, fg_color="transparent")
        frame_ImportCosine_std_csv.pack(fill='x', padx=10, pady=5)
        # Section for Standards CSV file input
        label_cosine_std_csv = ctk.CTkLabel(frame_ImportCosine_std_csv, text="Import Standards MS Spectra (CSV file):", font=('Arial', 14, 'bold'))
        label_cosine_std_csv.pack(side='left', padx=(0, 10))
        self.entry_cosine_std_csv_path = ctk.CTkEntry(frame_ImportCosine_std_csv, placeholder_text="Input Standards CSV file path", width=385)
        self.entry_cosine_std_csv_path.pack(side='left', padx=(0, 5))
        self.button_cosine_std_csv_browse = ctk.CTkButton(frame_ImportCosine_std_csv, text="Browse file",
                                            command=self.browse_std_csv, fg_color="#6598D3", hover_color="#5578B3")
        self.button_cosine_std_csv_browse.pack(padx=10, pady=5, anchor='e')

        # mass range of product ions
        frame_Cosine_para1 = ctk.CTkFrame(frame_inside_canvas, fg_color="transparent")
        frame_Cosine_para1.pack(fill='x', padx=10, pady=5)

        label_massrange_min  = ctk.CTkLabel(frame_Cosine_para1, text="・Mass Range of Product Ions (m/z):", font=('Arial', 14))
        label_massrange_min.pack(side='left', padx=(0, 10))
        self.entry_mass_range_min_Cosine = ctk.CTkEntry(frame_Cosine_para1, placeholder_text="min, e.g.200", width=203)
        self.entry_mass_range_min_Cosine.pack(side='left', padx=(0, 5))

        label_massrange_max  = ctk.CTkLabel(frame_Cosine_para1, text="to", font=('Arial', 14))
        label_massrange_max.pack(side='left', padx=(0, 5))

        # Option Menu for max_mass
        self.max_mass_option = ctk.CTkOptionMenu(frame_Cosine_para1, values=["Use Precursor Ion", "Specify Value"],
                                                fg_color="#F5F5F5",  # Background color
                                                button_color="#BEBEBE",  # Button color
                                                text_color="#333333",  # Text color
                                                width=160)
        self.max_mass_option.pack(side='left', padx=(0, 5))
        self.entry_mass_range_max_Cosine = ctk.CTkEntry(frame_Cosine_para1, placeholder_text="max, e.g.400", width=190)
        self.entry_mass_range_max_Cosine.pack(side='left', padx=(0, 5))

        # feature_formulas; mass range of product ions
        frame_Cosine_para2 = ctk.CTkFrame(frame_inside_canvas, fg_color="transparent")
        frame_Cosine_para2.pack(fill='x', padx=10, pady=5)
        label_IntensityNor_Cosine  = ctk.CTkLabel(frame_Cosine_para2, text="・Intensity Normalization:", font=('Arial', 14))
        label_IntensityNor_Cosine.pack(side='left', padx=(0, 5))
        self.entry_IntensityNor_Cosine = ctk.CTkEntry(frame_Cosine_para2, placeholder_text="e.g., 0.1", width=271)
        self.entry_IntensityNor_Cosine.pack(side='left', padx=(0, 5))

        label_Cosine  = ctk.CTkLabel(frame_Cosine_para2, text="・Cosine Score:", font=('Arial', 14))
        label_Cosine.pack(side='left', padx=(0, 5))
        self.entry_CosineScore = ctk.CTkEntry(frame_Cosine_para2, placeholder_text="e.g., 0.7", width=271)
        self.entry_CosineScore.pack(side='left', padx=(0, 5))

        # CTkFrame of Output
        frame_Output = ctk.CTkFrame(frame_inside_canvas, fg_color="transparent")
        frame_Output.pack(fill='x', padx=10, pady=5)
        # Section for Output file settings
        label_output = ctk.CTkLabel(frame_Output, text="Export Result Files:", font=('Arial', 14, 'bold'))
        label_output.pack(side='left', padx=(0, 10))
        self.entry_cosine_output_path = ctk.CTkEntry(frame_Output,
                                                        placeholder_text="Output directory path, e.g. C:/Users/xxx/", width=315)
        self.entry_cosine_output_path.pack(side='left', padx=(0, 5))
        self.entry_cosine_output_csv_filename = ctk.CTkEntry(frame_Output,
                                                    placeholder_text="Output CSV file name, e.g. result.csv", width=205)
        self.entry_cosine_output_csv_filename.pack(side='left', padx=(0, 5))
        self.button_cosine_calculate = ctk.CTkButton(frame_Output, text="Calculate",
                                                        command=self.calculate_cosine, fg_color="#6598D3", hover_color="#5578B3")
        self.button_cosine_calculate.pack(padx=10, pady=5, anchor='e')

    def browse_query_mgf(self):
        path = filedialog.askopenfilename(filetypes=[("MGF files", "*.mgf")])
        if path:
            self.entry_cosine_query_mgf_path.delete(0, tk.END)
            self.entry_cosine_query_mgf_path.insert(0, path)

    def browse_query_csv(self):
        path = filedialog.askopenfilename(filetypes=[("CSV files", "*.csv")])
        if path:
            self.entry_cosine_query_csv_path.delete(0, tk.END)
            self.entry_cosine_query_csv_path.insert(0, path)

    def browse_std_mgf(self):
        path = filedialog.askopenfilename(filetypes=[("MGF files", "*.mgf")])
        if path:
            self.entry_cosine_std_mgf_path.delete(0, tk.END)
            self.entry_cosine_std_mgf_path.insert(0, path)

    def browse_std_csv(self):
        path = filedialog.askopenfilename(filetypes=[("CSV files", "*.csv")])
        if path:
            self.entry_cosine_std_csv_path.delete(0, tk.END)
            self.entry_cosine_std_csv_path.insert(0, path)

    def calculate_cosine(self):
        try:
            output_dir = self.entry_cosine_output_path.get().strip()
            output_filename = self.entry_cosine_output_csv_filename.get().strip()
            full_output_path = os.path.join(output_dir, output_filename)

            if not os.path.exists(output_dir):
                os.makedirs(output_dir)

            self.process_cosine_score_calculate(
                self.entry_cosine_query_mgf_path.get(),
                self.entry_cosine_std_mgf_path.get(),
                self.entry_cosine_query_csv_path.get(),
                self.entry_cosine_std_csv_path.get(),
                full_output_path
            )
            messagebox.showinfo("Process Complete", "Cosine scores have been successfully calculated.")
        except Exception as e:
            messagebox.showerror("Error", f"Failed to save file: {str(e)}")
            print(f"Failed to save file: {e}")

    def filter_spectra_with_normalization(self, filename):

        try:
            min_mass = float(self.entry_mass_range_min_Cosine.get())
            norScore = float(self.entry_IntensityNor_Cosine.get())
            max_mass_option = self.max_mass_option.get()
            if max_mass_option == "Specify Value":
                max_mass = float(self.entry_mass_range_max_Cosine.get())
        except ValueError:
            messagebox.showerror("Input Error",
                                "Please enter valid numeric values")
            return {}

        filtered_spectra = []
        with mgf.read(filename) as spectra:
            for spectrum in spectra:
                mz_array = spectrum['m/z array']
                intensity_array = spectrum['intensity array']
                if max_mass_option == "Use Precursor Ion":
                    max_mass = spectrum['params']['pepmass'][0]
                ind = (mz_array >= min_mass) & (mz_array <= max_mass) & (intensity_array > norScore)
                if np.any(ind):
                    mz_array = mz_array[ind]
                    intensity_array = intensity_array[ind]
                    intensity_max = np.max(intensity_array)
                    intensity_min = np.min(intensity_array)
                    normalized_intensities = (intensity_array - intensity_min) / (intensity_max - intensity_min)
                    filtered_spectra.append({
                        'm/z array': mz_array,
                        'intensity array': normalized_intensities,
                        'params': spectrum['params']
                    })
        return filtered_spectra

    def convert_to_matchms_spectra(self, filtered_spectra):
        matchms_spectra = []
        for spec in filtered_spectra:
            spectrum = Spectrum(mz=spec['m/z array'],
                                intensities=spec['intensity array'],
                                metadata=spec['params'])
            matchms_spectra.append(spectrum)
        return matchms_spectra

    def process_cosine_score_calculate(self, query_path_mgf, std_path_mgf, query_path_csv, std_path_csv, full_output_path):

        try:
            CosineScore = float(self.entry_CosineScore.get())
        except ValueError:
            messagebox.showerror("Input Error",
                                "Please enter valid numeric values")
            return {}

        query_spectra_filtered = self.convert_to_matchms_spectra(self.filter_spectra_with_normalization(query_path_mgf))
        std_spectra_filtered = self.convert_to_matchms_spectra(self.filter_spectra_with_normalization(std_path_mgf))

        similarity_scores = calculate_scores(query_spectra_filtered, std_spectra_filtered,
                                            CosineGreedy(tolerance=0.005), is_symmetric=False)

        query_csv = pd.read_csv(query_path_csv)
        std_csv = pd.read_csv(std_path_csv)
        query_csv['row ID'] = query_csv['row ID'].astype(str)
        std_csv['row ID'] = std_csv['row ID'].astype(str)

        best_matches = {}
        cosine_scores = {}

        for score in similarity_scores:
            query_spectrum, std_spectrum, score_values = score
            cosine_score, matched_peaks = score_values
            if cosine_score > CosineScore:
                query_title = query_spectrum.metadata['title']
                std_title = std_spectrum.metadata['title']
                if std_title in std_csv['row ID'].values:
                    b_ring_structure = std_csv.loc[std_csv['row ID'] == std_title, 'B ring structure'].iloc[0]
                    best_matches[query_title] = b_ring_structure
                    cosine_scores[query_title] = cosine_score

        query_csv['cosine score'] = query_csv['row ID'].map(cosine_scores).round(2).fillna('value < 0.5')
        query_csv['B part'] = query_csv['row ID'].map(best_matches).fillna('others')

        insert_index = query_csv.columns.get_loc('C7 loss') + 1
        columns_to_insert = ['cosine score', 'B part']
        for i, col in enumerate(columns_to_insert):
            query_csv.insert(insert_index + i, col, query_csv.pop(col))

        query_csv.reset_index(inplace=True)
        query_csv.to_csv(full_output_path, index=False)


    '''annotation_tab'''
    def create_annotation_tab(self):
        annotation_tab = self.tab("7.Annotation")

        # Create a canvas to contain all the frames
        canvas = tk.Canvas(annotation_tab)
        canvas.pack(side="left", fill="both", expand=True)

        # Create a main frame inside the canvas
        main_frame = ctk.CTkFrame(canvas)
        canvas.create_window((0, 0), window=main_frame, anchor='nw')

        # Part I
        # CTkFrame of Section label
        frame_annotation_sectionI_label = ctk.CTkFrame(main_frame, fg_color="transparent")
        frame_annotation_sectionI_label.pack(fill='x', padx=10, pady=5)
        # Section label
        label_sectionI = ctk.CTkLabel(frame_annotation_sectionI_label, text="Annotation Part I",
                                    font=('Arial', 16, 'bold'))
        label_sectionI.pack(side='left', padx=(0, 10))

        # CTkFrame of Import CSV
        frame_annotationI_query_csv = ctk.CTkFrame(main_frame, fg_color="transparent")
        frame_annotationI_query_csv.pack(fill='x', padx=10, pady=5)
        # Section for Query CSV file input
        label_annotationI_query_csv = ctk.CTkLabel(frame_annotationI_query_csv,
                                                text="  Import Query MS Spectra (CSV file):",
                                                font=('Arial', 14, 'bold'))
        label_annotationI_query_csv.pack(side='left', padx=(0, 10))
        self.entry_annotationI_query_csv_path = ctk.CTkEntry(frame_annotationI_query_csv,
                                                            placeholder_text="Input query CSV file path", width=440)
        self.entry_annotationI_query_csv_path.pack(side='left', padx=(0, 5))
        self.button_annotationI_query_csv_browse = ctk.CTkButton(frame_annotationI_query_csv, text="Browse file",
                                                                command=self.browse_annotationPartI_query_csv, fg_color="#6598D3", hover_color="#5578B3")
        self.button_annotationI_query_csv_browse.pack(padx=10, pady=5, anchor='e')

        # CTkFrame of Import acylation CSV
        frame_annotationI_Acyl_csv = ctk.CTkFrame(main_frame, fg_color="transparent")
        frame_annotationI_Acyl_csv.pack(fill='x', padx=10, pady=5)
        # Section for Import acylation CSV
        label_annotationI_Acyl_csv = ctk.CTkLabel(frame_annotationI_Acyl_csv,
                                                text="  Import Simulated Acylation file (CSV file):",
                                                font=('Arial', 14, 'bold'))
        label_annotationI_Acyl_csv.pack(side='left', padx=(0, 10))
        self.entry_annotationI_Acyl_csv_path = ctk.CTkEntry(frame_annotationI_Acyl_csv,
                                                            placeholder_text="Input Simulated Acylation CSV file path",
                                                            width=400)
        self.entry_annotationI_Acyl_csv_path.pack(side='left', padx=(0, 5))
        self.button_annotationI_Acyl_csv_browse = ctk.CTkButton(frame_annotationI_Acyl_csv, text="Browse file",
                                                                command=self.browse_annotationAcylPartI_csv, fg_color="#6598D3", hover_color="#5578B3")
        self.button_annotationI_Acyl_csv_browse.pack(padx=10, pady=5, anchor='e')

        # CTkFrame of Output
        frame_annotationI_Output = ctk.CTkFrame(main_frame, fg_color="transparent")
        frame_annotationI_Output.pack(fill='x', padx=10, pady=5)
        # Section for Output file settings
        label_annotationI_output = ctk.CTkLabel(frame_annotationI_Output, text="  Export Result Files:",
                                                font=('Arial', 14, 'bold'))
        label_annotationI_output.pack(side='left', padx=(0, 10))
        self.entry_annotationI_output_path = ctk.CTkEntry(frame_annotationI_Output,
                                                        placeholder_text="Output directory path, e.g. C:/Users/xxx/",
                                                        width=320)
        self.entry_annotationI_output_path.pack(side='left', padx=(0, 5))
        self.entry_annotationI_output_csv_filename = ctk.CTkEntry(frame_annotationI_Output,
                                                                placeholder_text="Output CSV file name, e.g. result.csv",
                                                                width=228)
        self.entry_annotationI_output_csv_filename.pack(side='left', padx=(0, 5))
        self.button_annotationI_calculate = ctk.CTkButton(frame_annotationI_Output, text="Annotate Part I",
                                                        command=self.annotateI, fg_color="#6598D3", hover_color="#5578B3")
        self.button_annotationI_calculate.pack(padx=10, pady=5, anchor='e')

        # Part II
        # CTkFrame of Section label
        frame_annotation_sectionII_label = ctk.CTkFrame(main_frame, fg_color="transparent")
        frame_annotation_sectionII_label.pack(fill='x', padx=10, pady=20)
        # Section label
        label_sectionII = ctk.CTkLabel(frame_annotation_sectionII_label, text="Annotation Part II",
                                    font=('Arial', 16, 'bold'))
        label_sectionII.pack(side='left', padx=(0, 10))

        # CTkFrame of Import CSV
        frame_annotationII_query_csv = ctk.CTkFrame(main_frame, fg_color="transparent")
        frame_annotationII_query_csv.pack(fill='x', padx=10, pady=5)
        # Section for Query CSV file input
        label_annotationII_query_csv = ctk.CTkLabel(frame_annotationII_query_csv,
                                                    text="  Import Query MS Spectra (CSV file):",
                                                    font=('Arial', 14, 'bold'))
        label_annotationII_query_csv.pack(side='left', padx=(0, 10))
        self.entry_annotationII_query_csv_path = ctk.CTkEntry(frame_annotationII_query_csv,
                                                            placeholder_text="Input query CSV file path", width=440)
        self.entry_annotationII_query_csv_path.pack(side='left', padx=(0, 5))
        self.button_annotationII_query_csv_browse = ctk.CTkButton(frame_annotationII_query_csv, text="Browse file",
                                                                command=self.browse_annotationPartII_query_csv, fg_color="#6598D3", hover_color="#5578B3")
        self.button_annotationII_query_csv_browse.pack(padx=10, pady=5, anchor='e')

        # CTkFrame of Output
        frame_annotationII_Output = ctk.CTkFrame(main_frame, fg_color="transparent")
        frame_annotationII_Output.pack(fill='x', padx=10, pady=5)
        # Section for Output file settings
        label_annotationII_output = ctk.CTkLabel(frame_annotationII_Output, text="  Export Result Files:",
                                                font=('Arial', 14, 'bold'))
        label_annotationII_output.pack(side='left', padx=(0, 10))
        self.entry_annotationII_output_path = ctk.CTkEntry(frame_annotationII_Output,
                                                        placeholder_text="Output directory path, e.g. C:/Users/xxx/",
                                                        width=320)
        self.entry_annotationII_output_path.pack(side='left', padx=(0, 5))
        self.entry_annotationII_output_csv_filename = ctk.CTkEntry(frame_annotationII_Output,
                                                                placeholder_text="Output CSV file name, e.g. result.csv",
                                                                width=228)
        self.entry_annotationII_output_csv_filename.pack(side='left', padx=(0, 5))
        self.button_annotationII_calculate = ctk.CTkButton(frame_annotationII_Output, text="Annotate Part II",
                                                        command=self.annotateII, fg_color="#6598D3", hover_color="#5578B3")
        self.button_annotationII_calculate.pack(padx=10, pady=5, anchor='e')

    def browse_annotationPartI_query_csv(self):
        path = filedialog.askopenfilename(filetypes=[("CSV files", "*.csv")])
        if path:
            self.entry_annotationI_query_csv_path.delete(0, tk.END)
            self.entry_annotationI_query_csv_path.insert(0, path)

    def browse_annotationAcylPartI_csv(self):
        path = filedialog.askopenfilename(filetypes=[("CSV files", "*.csv")])
        if path:
            self.entry_annotationI_Acyl_csv_path.delete(0, tk.END)
            self.entry_annotationI_Acyl_csv_path.insert(0, path)

    def browse_annotationPartII_query_csv(self):
        path = filedialog.askopenfilename(filetypes=[("CSV files", "*.csv")])
        if path:
            self.entry_annotationII_query_csv_path.delete(0, tk.END)
            self.entry_annotationII_query_csv_path.insert(0, path)

    def annotateI(self):
        try:
            output_dir = self.entry_annotationI_output_path.get().strip()
            output_filename = self.entry_annotationI_output_csv_filename.get().strip()
            full_output_path = os.path.join(output_dir, output_filename)

            if not os.path.exists(output_dir):
                os.makedirs(output_dir)

            df = pd.read_csv(self.entry_annotationI_query_csv_path.get())
            acyl_df = pd.read_csv(self.entry_annotationI_Acyl_csv_path.get())

            def count_element(molecular_formula, element):
                pattern = re.compile(rf'{element}(\d*)')
                match = pattern.search(molecular_formula)
                if match:
                    return int(match.group(1) or 1)
                return 0

            def determine_a_parts(row):
                if row['type'] == 'D':
                    o_count = count_element(row['C17 Ion'], 'O')
                    h_count = count_element(row['C17 Ion'], 'H')
                    o_count_c18 = count_element(row.get('C18 Ion', ''), 'O')
                    if o_count == 2:
                        if h_count == 16:
                            return '1,2-en-3-one'
                        elif h_count == 20:
                            return '1,2-dihydro-3-ol'
                        elif h_count == 18:
                            if o_count_c18 == 1:
                                return '1,2-en-3-ol'
                            elif o_count_c18 == 2:
                                return '1,2-dihydro-3-one'
                    elif o_count == 3:
                        if h_count == 16:
                            return '1,2-en-3-one'
                        elif h_count == 20:
                            return '1,2-dihydro-3-ol'
                        elif h_count == 18:
                            if o_count_c18 == 1:
                                return '1,2-en-3-ol'
                            elif o_count_c18 == 2:
                                return '1,2-dihydro-3-one'
                elif row['type'] == 'MD':
                    if row['C15 Ion'] != 'none' and row.get('C3 loss', '') != 'none':
                        return 'bicyclo [2.2.1] heptane ring'
                    elif row['C15 Ion'] == 'none':
                        return '1-alkyl-3-ol'
                return 'others'

            df['A part'] = df.apply(determine_a_parts, axis=1)

            def determine_b_parts(row):
                if row['type'] == 'D':
                    return row['B part']
                elif row['type'] == 'MD':
                    c_count_c27 = count_element(row.get('C27 Ion', ''), 'C')
                    if c_count_c27 == 27 or row.get('C3 loss', '') != 'none':
                        return '6,7-epoxy'
                    elif row.get('C27 Ion', '') == 'none' or row.get('C3 loss', '') == 'none':
                        return 'others'

            df['B part'] = df.apply(determine_b_parts, axis=1)

            def determine_c12(row):
                if row['type'] == 'D':
                    o_count = count_element(row['C17 Ion'], 'O')
                    if o_count == 2:
                        return 'no C-12'
                    elif o_count == 3:
                        return 'with C-12'
                return 'others'

            df['C-12'] = df.apply(determine_c12, axis=1)

            def determine_m_parts(row):
                if row['type'] == 'MD':
                    h_count = count_element(row.get('C10 loss', ''), 'H')
                    if h_count == 18:
                        return 'C10 ring with no substituents'
                    elif h_count == 16:
                        return 'C10 ring with one substituents'
                    elif h_count == 14:
                        return 'C10 ring with two substituents'
                    elif h_count == 12:
                        return 'C10 ring with three substituents'
                    elif h_count == 10:
                        return 'C10 ring with four substituents'
                    c_count = count_element(row.get('C14 loss', ''), 'C')
                    if c_count == 14:
                        return 'C14 ring with an olenfinic'
                return 'others'

            df['M part'] = df.apply(determine_m_parts, axis=1)

            def determine_acyl_parts(row):
                acyl_parts = []
                if row['C2 loss'] == 'C2H4O2':
                    acyl_parts.append('acyl')

                if row['Acyl Ion'] != 'none':
                    acyl_ions = row['Acyl Ion'].split(',')
                    acyl_formulas = acyl_df['Acyl Formula'].tolist()
                    identifications = acyl_df['Identification'].tolist()

                    for ion in acyl_ions:
                        matches = [identifications[i] for i, formula in enumerate(acyl_formulas) if formula.strip() == ion.strip()]
                        if matches:
                            acyl_parts.extend(matches)

                if not acyl_parts:
                    return 'others'
                return ','.join(acyl_parts)

            df['Acyl part'] = df.apply(determine_acyl_parts, axis=1)

            def determine_acyl_numb(row):
                acyl_numb = []
                if row['C2 loss'] == 'C2H4O2':
                    acyl_numb.append('R1')

                if row['Acyl Ion'] != 'none':
                    acyl_ions = row['Acyl Ion'].split(',')
                    acyl_formulas = acyl_df['Acyl Formula'].tolist()
                    r_parts = acyl_df['R parts'].tolist()

                    for ion in acyl_ions:
                        matches = [r_parts[i] for i, formula in enumerate(acyl_formulas) if formula.strip() == ion.strip()]
                        if matches:
                            acyl_numb.extend(matches)

                if not acyl_numb:
                    return 'others'
                return ','.join(acyl_numb)

            df['Acyl Number'] = df.apply(determine_acyl_numb, axis=1)

            def determine_acyl_simile(row):
                acyl_simile = []
                if row['C2 loss'] == 'C2H4O2':
                    acyl_simile.append('CC(=O)*')

                if row['Acyl Ion'] != 'none':
                    acyl_ions = row['Acyl Ion'].split(',')
                    acyl_formulas = acyl_df['Acyl Formula'].tolist()
                    acyl_SMILES = acyl_df['Acyl SMILES'].tolist()

                    for ion in acyl_ions:
                        matches = [acyl_SMILES[i] for i, formula in enumerate(acyl_formulas) if formula.strip() == ion.strip()]
                        if matches:
                            acyl_simile.extend(matches)

                if not acyl_simile:
                    return 'others'
                return ','.join(acyl_simile)

            df['Acyl SMILES'] = df.apply(determine_acyl_simile, axis=1)

            df.to_csv(full_output_path, index=False)

            messagebox.showinfo("Process Complete", "Annotation has been successfully conducted.")
        except Exception as e:
            messagebox.showerror("Error", f"Failed to save file: {str(e)}")
            print(f"Failed to save file: {e}")

    def annotateII(self):
        try:
            output_dir = self.entry_annotationII_output_path.get().strip()
            output_filename = self.entry_annotationII_output_csv_filename.get().strip()
            full_output_path = os.path.join(output_dir, output_filename)

            if not os.path.exists(output_dir):
                os.makedirs(output_dir)

            df = pd.read_csv(self.entry_annotationII_query_csv_path.get())

            data_a = {
                'A number': ['A1', 'A2', 'A3', 'A4', 'A5', 'A6', 'A7', 'A8', 'A9', 'A10', 'A11'],
                'A Struc': [
                    '1,2-en-3-one', '1,2-dihydro-3-one', '1,2-en-3-ol', '1,2-dihydro-3-ol', '1,10-en-3-one',
                    '1-alkyl-3-one', '1-alkyl-2-ol-3-one', '1-alkyl-3-ol', '1-alkyl-2-ol-3-ol',
                    '3,4-seco ring', 'bicyclo [2.2.1] heptane ring'
                ]
            }
            data_b = {
                'B number': ['B1', 'B2', 'B3', 'B4', 'B5', 'B6', 'B7', 'B8'],
                'B Struc': [
                    '6,7-epoxy', '6,7-diol', '6,7-ene',
                    '5-dehydro-6,7-epoxy', '5-dehydro-6,7-diol', '5-dehydro-6,7-ene',
                    '4,7-epoxy', '4,6-epoxy'
                ]
            }
            data_c = {
                'C number': ['C1', 'C2', 'C3', 'C4', 'C5', 'C6'],
                'C Struc': [
                    '9, 13, 14-orthoester', '9, 13, 14-orthoester', '9, 13, 14-triol', '9, 13, 14-triol',
                    '9, 13, 14-orthoester and 18-ol', '9, 13, 14-orthoester and 18-ol'
                ],
                'C-12Struc': [
                    'no C-12', 'with C-12', 'no C-12', 'with C-12', 'no C-12', 'with C-12'
                ]
            }
            data_m = {
                'M number': ['M1', 'M2', 'M3', 'M4', 'M5', 'M6', 'M7', 'M8', 'M9'],
                'M Struc': [
                    'C9 ring with no substituents', 'C10 ring with no substituents',
                    'C10 ring with one substituents', 'C10 ring with two substituents',
                    'C10 ring with three substituents', 'C10 ring with four substituents',
                    'C14 ring with an olenfinic',
                    'C16 ring with an olenfinic', 'C16 ring with one substituents'
                ]
            }

            df_a = pd.DataFrame(data_a)
            df_b = pd.DataFrame(data_b)
            df_c = pd.DataFrame(data_c)
            df_m = pd.DataFrame(data_m)

            # A number
            def map_a_number(row):
                match = df_a[df_a['A Struc'] == row['A part']]
                if not match.empty:
                    return match['A number'].values[0]
                return 'others'

            df['A number'] = df.apply(map_a_number, axis=1)

            # B number
            def map_b_number(row):
                match = df_b[df_b['B Struc'] == row['B part']]
                if not match.empty:
                    return match['B number'].values[0]
                return 'others'

            df['B number'] = df.apply(map_b_number, axis=1)

            # C number
            def map_c_number(row):
                if row['type'] == 'D':
                    match = df_c[(df_c['C Struc'] == row['C part']) & (df_c['C-12Struc'] == row['C-12'])]
                    if not match.empty:
                        return match['C number'].values[0]
                elif row['type'] == 'MD':
                    match = df_c[df_c['C Struc'] == row['C part']]
                    if not match.empty:
                        return match['C number'].values[0]
                return 'others'

            df['C number'] = df.apply(map_c_number, axis=1)

            # M number
            def map_m_number(row):
                match = df_m[df_m['M Struc'] == row['M part']]
                if not match.empty:
                    return match['M number'].values[0]
                return 'others'

            df['M number'] = df.apply(map_m_number, axis=1)

            # Insert the new columns after 'Acyl part'
            insert_index = df.columns.get_loc('Acyl part') + 1
            for col in ['A number', 'B number', 'C number', 'M number']:
                cols = df.columns.tolist()
                cols.insert(insert_index, cols.pop(cols.index(col)))
                df = df[cols]
                insert_index += 1

            # Save the updated DataFrame to a new CSV file
            df.to_csv(full_output_path, index=False)

            messagebox.showinfo("Process Complete", "Annotation has been successfully conducted.")
        except Exception as e:
            messagebox.showerror("Error", f"Failed to save file: {str(e)}")
            print(f"Failed to save file: {e}")


    '''visualization_tab'''
    def create_visualization_tab(self):
        visualization_tab = self.tab("8.Visualization")

        frame_visualization_query_csv = ctk.CTkFrame(visualization_tab)
        frame_visualization_query_csv.pack(fill='x', padx=10, pady=5)
        label_visualization_query_csv = ctk.CTkLabel(frame_visualization_query_csv, text="Import Query Metadata (CSV file):",
                                                font=('Arial', 14, 'bold'))
        label_visualization_query_csv.pack(side='left', padx=(0, 10))
        self.entry_visualization_query_csv_path = ctk.CTkEntry(frame_visualization_query_csv,
                                                            placeholder_text="Input query Metadata file path",
                                                            width=570)
        self.entry_visualization_query_csv_path.pack(side='left', padx=(0, 5))
        self.button_visualization_query_csv_browse = ctk.CTkButton(frame_visualization_query_csv, text="Browse file",
                                                                command=self.browse_visualization_query_csv, fg_color="#6598D3", hover_color="#5578B3")
        self.button_visualization_query_csv_browse.pack(padx=10, pady=5, anchor='e')

        frame_visualization_SkeletonDB_csv = ctk.CTkFrame(visualization_tab)
        frame_visualization_SkeletonDB_csv.pack(fill='x', padx=10, pady=5)
        label_visualization_SkeletonDB_csv = ctk.CTkLabel(frame_visualization_SkeletonDB_csv,
                                                    text="Import Skeleton Database (CSV file):",
                                                    font=('Arial', 14, 'bold'))
        label_visualization_SkeletonDB_csv.pack(side='left', padx=(0, 10))
        self.entry_visualization_SkeletonDB_csv_path = ctk.CTkEntry(frame_visualization_SkeletonDB_csv,
                                                                placeholder_text="Input db CSV file path", width=550)
        self.entry_visualization_SkeletonDB_csv_path.pack(side='left', padx=(0, 5))
        self.button_visualization_SkeletonDB_csv_browse = ctk.CTkButton(frame_visualization_SkeletonDB_csv,
                                                                    text="Browse file",
                                                                    command=self.browse_visualization_SkeletonDB_csv, fg_color="#6598D3", hover_color="#5578B3")
        self.button_visualization_SkeletonDB_csv_browse.pack(padx=10, pady=5, anchor='e')

        frame_Output = ctk.CTkFrame(visualization_tab)
        frame_Output.pack(fill='x', padx=10, pady=5)
        label_output = ctk.CTkLabel(frame_Output, text="Export Result Files:", font=('Arial', 14, 'bold'))
        label_output.pack(side='left', padx=(0, 10))
        self.entry_visualization_output_path = ctk.CTkEntry(frame_Output,
                                                        placeholder_text="Output directory path, e.g. C:/Users/xxx/",
                                                        width=415)
        self.entry_visualization_output_path.pack(side='left', padx=(0, 5))
        self.entry_visualization_output_csv_filename = ctk.CTkEntry(frame_Output,
                                                                placeholder_text="Output CSV file name, e.g. result.csv",
                                                                width=245)
        self.entry_visualization_output_csv_filename.pack(side='left', padx=(0, 5))
        self.button_visualization_calculate = ctk.CTkButton(frame_Output, text="Calculate",
                                                        command=self.calculate_visualization, fg_color="#6598D3", hover_color="#5578B3")
        self.button_visualization_calculate.pack(padx=10, pady=5, anchor='e')

        # Add image display panels on the right
        self.structure_display_frame = ctk.CTkFrame(visualization_tab, width=400)
        self.structure_display_frame.pack(side='right', fill='y', expand=False, padx=10, pady=5)

        self.core_structure_label = ctk.CTkLabel(self.structure_display_frame, text="Core Structure",
                                                font=('Arial', 14, 'bold'))
        self.core_structure_label.pack(pady=(10, 5))
        self.core_structure_canvas = tk.Canvas(self.structure_display_frame, width=380, height=300)
        self.core_structure_canvas.pack()

        self.substituents_structure_label = ctk.CTkLabel(self.structure_display_frame, text="Substituents Structure",
                                                font=('Arial', 14, 'bold'))
        self.substituents_structure_label.pack(pady=(10, 5))
        self.substituents_structure_canvas = tk.Canvas(self.structure_display_frame, width=380,
                                                    height=300)
        self.substituents_structure_canvas.pack()

        # Treeview for displaying CSV data
        self.treeview_frame = ctk.CTkFrame(visualization_tab)
        self.treeview_frame.pack(fill='both', expand=True, padx=10, pady=5)

        self.tree = ttk.Treeview(self.treeview_frame)
        self.tree.pack(fill='both', expand=True, side='left')
        self.tree.bind("<<TreeviewSelect>>", self.update_structure_displays)

        self.tree_scroll_y = tk.Scrollbar(self.treeview_frame, orient='vertical', command=self.tree.yview)
        self.tree_scroll_y.pack(side='right', fill='y')
        self.tree.configure(yscrollcommand=self.tree_scroll_y.set)

        self.tree_scroll_x = tk.Scrollbar(self.treeview_frame, orient='horizontal', command=self.tree.xview)
        self.tree_scroll_x.pack(side='bottom', fill='x')
        self.tree.configure(xscrollcommand=self.tree_scroll_x.set)

        self.substituents_images = []

        self.scale_frame = ctk.CTkFrame(visualization_tab)
        self.scale_frame.pack(fill='x', padx=10, pady=5)
        self.zoom_label = ctk.CTkLabel(self.scale_frame, text="Zoom:")
        self.zoom_label.pack(side='left', padx=(10, 5))
        self.zoom_scale = ctk.CTkSlider(self.scale_frame, from_=50, to=200, command=self.zoom_tree)
        self.zoom_scale.set(100)
        self.zoom_scale.pack(side='left', fill='x', expand=True, padx=(5, 10))
        self.zoom_percentage_label = ctk.CTkLabel(self.scale_frame, text="100%")
        self.zoom_percentage_label.pack(side='left', padx=(5, 10))

        self.default_font_size = 10

    def smiles_to_image(self, smiles, canvas, x_offset=10, label=None):
        try:
            mol = Chem.MolFromSmiles(smiles)
            if mol is None:
                return
            img = Draw.MolToImage(mol, size=(300, 200))
            image = ImageTk.PhotoImage(image=img)
            canvas.create_image(x_offset, 10, anchor='nw', image=image)
            self.substituents_images.append(image)
            if label:
                canvas.create_text(x_offset + 150, 180, text=label, font=('Arial', 12, 'bold'))
        except Exception as e:
            print(f"Error converting SMILES to image: {smiles}, error: {e}")

    def update_structure_displays(self, event):
        selected_item = self.tree.item(self.tree.focus())
        if not selected_item or 'values' not in selected_item:
            return

        skeleton_smiles_index = 36
        acyl_smiles_index = 37

        values = selected_item['values']
        if len(values) <= max(skeleton_smiles_index, acyl_smiles_index):
            return

        core_smiles = values[skeleton_smiles_index]
        substituent_smiles = values[acyl_smiles_index].split(',')

        self.core_structure_canvas.delete("all")
        if core_smiles != 'others':
            self.smiles_to_image(core_smiles, self.core_structure_canvas)

        self.substituents_structure_canvas.delete("all")
        x_offset = 10
        for i, smiles in enumerate(substituent_smiles):
            if smiles != 'others':
                self.smiles_to_image(smiles, self.substituents_structure_canvas, x_offset, label=f'R{i + 1}')
                x_offset += 150

    def browse_visualization_query_csv(self):
        path = filedialog.askopenfilename(filetypes=[("CSV files", "*.csv")])
        if path:
            self.entry_visualization_query_csv_path.delete(0, tk.END)
            self.entry_visualization_query_csv_path.insert(0, path)

    def browse_visualization_SkeletonDB_csv(self):
        path = filedialog.askopenfilename(filetypes=[("CSV files", "*.csv")])
        if path:
            self.entry_visualization_SkeletonDB_csv_path.delete(0, tk.END)
            self.entry_visualization_SkeletonDB_csv_path.insert(0, path)

    def calculate_visualization(self):
        try:
            query_csv_path = self.entry_visualization_query_csv_path.get().strip()
            skeleton_db_csv_path = self.entry_visualization_SkeletonDB_csv_path.get().strip()
            output_dir = self.entry_visualization_output_path.get().strip()
            output_filename = self.entry_visualization_output_csv_filename.get().strip()
            full_output_path = os.path.join(output_dir, output_filename)

            if not os.path.exists(output_dir):
                os.makedirs(output_dir)

            self.process_visualization(
                query_csv_path,
                skeleton_db_csv_path,
                full_output_path
            )
            self.display_csv_in_treeview(full_output_path)
            messagebox.showinfo("Process Complete", "Visualization process has been successfully completed.")
        except Exception as e:
            messagebox.showerror("Error", f"Failed to complete visualization: {str(e)}")
            print(f"Failed to complete visualization: {e}")

    def process_visualization(self, input_path_csv, skeleton_db_path_csv, output_path_csv):
        df = pd.read_csv(input_path_csv)
        skeleton_db = pd.read_csv(skeleton_db_path_csv)

        def find_skeleton_smiles(row, db):
            match = pd.DataFrame()

            if row['type'] == 'D':
                match = db[
                    (db['A parts'] == row['A number']) &
                    (db['B parts'] == row['B number']) &
                    (db['C parts'] == row['C number'])
                    ]
            elif row['type'] == 'MD':
                match = db[
                    (db['A parts'] == row['A number']) &
                    (db['B parts'] == row['B number']) &
                    (db['C parts'] == row['C number']) &
                    (db['M parts'] == row['M number'])
                    ]

            if not match.empty:
                return match['New SMILES'].values[0]
            return 'others'

        df['Skeleton SMILES'] = df.apply(find_skeleton_smiles, axis=1, db=skeleton_db)

        insert_index = df.columns.get_loc('Acyl Number') + 1

        cols = df.columns.tolist()
        cols.insert(insert_index, cols.pop(cols.index('Skeleton SMILES')))
        df = df[cols]

        df.to_csv(output_path_csv, index=False)

    def display_csv_in_treeview(self, csv_path):
        df = pd.read_csv(csv_path)
        self.tree.delete(*self.tree.get_children())
        self.tree["column"] = list(df.columns)
        self.tree["show"] = "headings"
        for column in self.tree["columns"]:
            self.tree.heading(column, text=column)

        df_rows = df.to_numpy().tolist()
        for row in df_rows:
            self.tree.insert("", "end", values=row)

    def zoom_tree(self, value):
        zoom_factor = int(value)
        self.zoom_percentage_label.configure(text=f"{zoom_factor}%")
        font_size = int(self.default_font_size * (zoom_factor / 100))
        style = ttk.Style()
        style.configure("Treeview", font=("TkDefaultFont", font_size))
        self.tree.tag_configure('font', font=("TkDefaultFont", font_size))


class MainApp(ctk.CTk):
    def __init__(self):
        super().__init__()

        self.title("CSNPs Structure Annotation Tool")
        self.geometry("1050x900")

        self.tab_view = MyTabView(master=self)
        self.tab_view.pack(padx=15, pady=15, fill='both', expand=True)

app = MainApp()
app.mainloop()