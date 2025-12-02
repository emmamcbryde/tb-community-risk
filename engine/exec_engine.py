# -*- coding: utf-8 -*-
"""
Spyder Editor.

This is a temporary script file.
"""
import os
print("Working directory set to:", os.getcwd())


from engine.params import extract_core_parameters, summarize_core_parameters

params = extract_core_parameters("data/parameters.csv")
summary = summarize_core_parameters(params)


print(summary)
