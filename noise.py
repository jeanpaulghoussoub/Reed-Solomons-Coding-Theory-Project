import random
from math import sqrt, exp, pi, sin, inf, log
import numpy.random


def make_random_noise(encoded_message, proba_error):     # probability of error
    out_encoded_message = list(encoded_message)
    for i in range(0, len(encoded_message)):
        if random.random() < proba_error:
            out_encoded_message[i] = make_error(encoded_message[i], 0, len(encoded_message) - 1)
    
    num_error = 0
    for i in range(0, len(encoded_message)):
        if encoded_message[i] != out_encoded_message[i]: num_error += 1
    return out_encoded_message, num_error


def make_error(value, min, max):
    # make random another value
    new_value = random.randint(min, max)
    while (new_value == value):
        new_value = random.randint(min, max)
    return new_value
