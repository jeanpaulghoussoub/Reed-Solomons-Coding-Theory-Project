from galois_field import *

class ReedSolomonError(Exception):
    pass

def make_field(prim, c_exp):
    global galois_field_log, galois_field_exp, field_char
    galois_field_log, galois_field_exp, field_char = init_tables(prim, c_exp)

def rs_find_error_evaluator(synd, err_loc, nsym):
    #Compute the error evaluator polynomial Omega
    _, remainder = galois_field_poly_div( galois_field_poly_mul(synd, err_loc), ([1] + [0]*(nsym + 1)) ) # Omega(x) = [ Synd(x) * Error_loc(x) ] mod x^(n-k+1)
    return remainder

def rs_calc_syndromes(msg, nsym, fcr):
    #Computes the syndromes polynomial
    synd = [0] * nsym
    for i in range(0, nsym):
        synd[i] = galois_field_poly_eval(msg, galois_field_pow(2,i+fcr))
    return [0] + synd # tu thap den cao

def rs_find_error_locator(synd, nsym):
    # compute error locator polynomial
    err_loc = [1]
    auxi_poly = [1] # auxiliary polynomial B(x)
    L = 0 # length of LFSR (Linear Feedback Shift Register)

    for r in range(1, nsym+1):
        # Compute delta
        delta = synd[r]
        for j in  range(1, len(err_loc)):
            delta = galois_field_sub( delta, galois_field_mul(err_loc[-(j+1)], synd[r-j]))
        
        auxi_poly += [0]

        if delta != 0:
            old_err_loc = err_loc
            err_loc = galois_field_poly_add(err_loc, galois_field_poly_scale(auxi_poly, delta))
            if (2*L <= r-1):
                L = r - L
                auxi_poly = galois_field_poly_scale(old_err_loc, galois_field_inverse(delta))
    errs = len(err_loc) - 1
    if errs * 2 > nsym:
        raise ReedSolomonError("Too many errors when compute error locator polynomial!")
    return err_loc # in tu bac cao xuong bac thap

def rs_find_errors(err_loc, nmess): # nmess is len(msg_in)
    # Find the roots of error polynomial by Chien's search.
    err_pos = [] # error postions
    for i in range(nmess):
        if galois_field_poly_eval(err_loc, galois_field_pow(2, i)) == 0: # if anpha^i is a root -->  err_pos = nmess - 1 - i
            err_pos.append(nmess - 1 - i)
    # Check the number of errors positions
    if len(err_pos) != len(err_loc) - 1:
        raise ReedSolomonError("Too many errors found by Chien Search!")
    return err_pos


def rs_correct_error(msg_in, synd, err_loc, err_pos, fcr):
    # computes the values (error magnitude) by Forney algorithm to correct message.

    # calculate errors evaluator polynomial Omega
    err_eval_poly = rs_find_error_evaluator(synd[::-1], err_loc, len(err_loc)-1)

    # get the error location polynomial X from the error positions in err_pos
    X = []
    roots = [len(msg_in) - 1 - pos for pos in err_pos]
    for i in roots:
        X.append(galois_field_pow(2, i))

    # Forney algorithm: compute the magnitudes
    E = [0] * (len(msg_in))
    for i, Xi in enumerate(X):
        Xi_inv = galois_field_inverse(Xi)

        # compute the denominator of the Forney's formula
        denominator = 1
        for j in range(len(X)):
            if j != i:
                factor_j = galois_field_sub(1, galois_field_mul(Xi_inv, X[j]))
                denominator = galois_field_mul(denominator, factor_j)
        # Check denominator == 0 ?
        if denominator == 0:
            raise ReedSolomonError("Could not find error magnitude!")

        # Compute the magnitude by the Forney's formula
        magnitude = galois_field_poly_eval(err_eval_poly, Xi_inv)     # Omega at x = Xi inverse
        magnitude = galois_field_mul(galois_field_pow(Xi, 1-fcr), magnitude)    
        magnitude = galois_field_div(magnitude, denominator) # magnitude value of the error, calculated by the Forney algorithm
        E[err_pos[i]] = magnitude
    msg_in = galois_field_poly_add(msg_in, E)
    return msg_in

def rs_correct_msg(msg_in, nsym, fcr):
    # decoding function

    if len(msg_in) > field_char: # message is too big
        raise ValueError("Message is too long!")

    msg_out = list(msg_in)

    # prepare the syndrome polynomial
    synd = rs_calc_syndromes(msg_out, nsym, fcr)
    if max(synd) == 0:
        return msg_out[:-nsym], msg_out[-nsym:]

    # compute the error locator polynomial using Berlekamp-Massey
    err_loc = rs_find_error_locator(synd, nsym)
    # find error postions in message using Chien search
    err_pos = rs_find_errors(err_loc[::-1] , len(msg_out))
    if err_pos is None:
        raise ReedSolomonError("Could not locate error!")

    # compute error magnitude polynomials, then correct errors
    msg_out = rs_correct_error(msg_out, synd, err_loc, err_pos, fcr)
    # check if successfully decoded message
    synd = rs_calc_syndromes(msg_out, nsym, fcr)
    if max(synd) > 0:
        raise ReedSolomonError("Could not correct message!")
    return msg_out[:-nsym], msg_out[-nsym:]