from read_csv import create_genes_values_dictionary
from math import sqrt
import numpy as np

NUM_OF_DISEASE = 3

def analyze(filename):
  genes_t_stats = {}
  genes_values = create_genes_values_dictionary(filename)
  for key, value in genes_values.items():
    genes_t_stats[key] = calculate_t_stat(value)

  print(genes_t_stats)

def calculate_t_stat(values):
  disease_part = values[0:NUM_OF_DISEASE]
  normal_part = values[NUM_OF_DISEASE:len(values)]
  disease_values = [num for num in disease_part if num != -99]
  normal_values = [num for num in normal_part if num != -99]

  return t_stat(disease_values, normal_values)


def t_stat(values1, values2):
  np_values1 = np.array(values1)
  np_values2 = np.array(values2)
  mean1 = np_values1.mean()
  mean2 = np_values2.mean()
  std1 = np.std(np_values1, ddof=1, dtype=np.float32)
  std2 = np.std(np_values2, ddof=1, dtype=np.float32)

  numerator = mean1 - mean2
  denominator = sqrt((std1 ** 2 / len(np_values1)) + (std2 ** 2 / len(np_values2)))

  return numerator / denominator

analyze("working_sample.csv")
