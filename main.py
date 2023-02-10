import test as ts
import encoder, decoder
import galois_field
from datetime import datetime

def main():
    # function for set code. For example, set_rs_code(3,2) to set RS(7,3)
    def set_rs_code(m, t):
        global n, k
        n = 2**m - 1
        k = n - 2*t
        print("RS(%d, %d)" % (n, k))
        exp_field = galois_field.galois_field_get_suitable_field(n)
        encoder.make_field(galois_field.galois_field_get_prim_of_field(exp_field), exp_field)
        decoder.make_field(galois_field.galois_field_get_prim_of_field(exp_field), exp_field)

    # Configuration of the parameters and input message
    b = 0
    m = 5
    t = 5

    set_rs_code(m, t)

    def test():
	    ts.make_one_test(n, k, b)

    test()


if __name__ == "__main__":
    main()

