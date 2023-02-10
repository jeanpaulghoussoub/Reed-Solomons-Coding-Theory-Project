import random
import noise
import decoder
import encoder
from math import sqrt, log10
from datetime import datetime

def generate_random_message(length, max_value):
    message = [0] * length
    for i in range(len(message)):
        message[i] = random.randint(1, max_value)
    return message

def make_one_test(n, k, b):
    message = generate_random_message(k, n)
    print("Sent random message: %s" % message)
    ecc_number = n - k
    #ecc_number = int(input("Enter number error correcting codes: "))
    # n = k + ecc_number # the remaining n-k symbols will be the ECC code (more is better)

    # Encoding the input message
    encoded_message = encoder.rs_encode_message(message, ecc_number, b)
    print("After encoding:		%s" % encoded_message)

    received_message, error_number = noise.make_random_noise(encoded_message, proba_error = 0.05)
    print("After making noise:	%s" % received_message)
    print("Number of error:	%d" % error_number)

    # Decoding/repairing the corrupted message
    corrected_message, corrected_ecc = decoder.rs_correct_message(received_message, ecc_number, b)
    print("Repaired:		%s" % (corrected_message+corrected_ecc))
    print("Restore message: %s" % corrected_message[:k])
