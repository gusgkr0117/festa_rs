from sympy.ntheory import factorint

d1 = 3**6 * 19**2 * 29**2 * 37**2 * 83**2 * 139**2 * 167**2 * 251**2 * \
    419**2 * 421**2 * 701**2 * 839**2 * 1009**2 * 1259**2 * 3061**2 * 3779**2
d2 = 5**4 * 7**3 * 11**2 * 13**2 * 17**2 * 41**2 * 43**2 * 71**2 * 89**2 * \
    127**2 * 211**2 * 281**2 * 503**2 * 631**2 * \
    2309**2 * 2521**2 * 2647**2 * 2729**2
dA = 59**2 * 3023**2 * 3359**2 * 4409**2 * 5039**2 * 6299**2 * \
    6719**2 * 9181**2 * 19531**2 * 22679**2 * 41161**2
dA1 = 59**2 * 6299**2 * 6719**2 * 9181**2
dA2 = 3023**2 * 3359**2 * 4409**2 * 5039**2 * 19531**2 * 22679**2 * 41161**2
p_array = [0xFFFFFFFFFFFFFFFF, 0xFFFFFFFFFFFFFFFF, 0xFFFFFFFFFFFFFFFF, 0xFFFFFFFFFFFFFFFF, 0xFFFFFFFFFFFFFFFF, 0xFFFFFFFFFFFFFFFF, 0xFFFFFFFFFFFFFFFF, 0xFFFFFFFFFFFFFFFF, 0xFFFFFFFFFFFFFFFF, 0x43FFFFFFFFFFFFFF,
           0x9AC3245C4D7BE6B3, 0x21D7DCD797059B7B, 0x8F19A73E323F6569, 0x841FED4773CFDB16, 0x02979D50DD13D09A, 0x01712922BAF59934, 0xBD1C756E54F72C15, 0xF6B3CF47C54370FE, 0xCEC87BD4C1480F2B, 0x11CF13E54B11406F, 0x000000000000176C]
Px = 0xb22cf14c3787cd20b3e1eada2117d2996a50298e79660009eedfa608aa257d9e8071259d6cbbc4ff5f9f871351240711d39cf94e4c18f3cc5a3fa905410682435f3201f334e7b86db6157567a513de96e94d04f4c29371ab1930e86c1850765979c7a4d84690dd347ccc91cc769ab73da65c8d0aac8160ac9ff398044d58e03b86de0e3b63074bec0ae9ba18c569d14412c7c30cff9c9295be9b2dafeece8954a26
Py = 0x522f3c01b38a431d9b5090140dd8b71cd578046c1fdd942f3f62e5d73b8e567a03cfd89660855a08bda905f8561a4f79aca23d08570ddcae0b5181a55aa14f763e2af6efdb934f17ec4690de842938d6d487c22fd699b440b459a1d9825932bec3f1eadd94c491708a5cccfb6b687006a4426df07e9edb6820930e54c3ed0e8784de5fbd874e2e731e233169daa5116a1beda4b2092e262bd3ba926641d45b8a118
Qx = 0x11238ff030c55a5c6dd642567069e2ca09eddfda10165cc7a20960022e3792fa436bcc5cf225a610eb25a54ba2f12be6825a194e42014a0e15f2503c9cce4755c2fcb3455f058cd7e1e1653751c72e982b4431ccd53267606a86cdcaa700a0eca73209516ea0d30f9713b46ee18b995b528cf58bece0a9404f35f9ec70ec0214422eec42c1d4b79fe4e0ca9f07d54a0cd6099ce149fd474bb80154fcdf4057795a20
Qy = 0x138fdc3012d1e48dd9048717fb8aee0c371ee0d4298f428d2ed71edd64cae0736797a611bad0d873b395a49d4b643d8b6d87ae5d9f859cd30a114a746a5a1fd93a438fa82803240603feb34dc5908c4abe6e1e3d6b1970db5f9d1995b5152d368056c6e44780fa7ee4a716239137be42d449069d2e0d2d3c5bf05fc9efaa795d975babf26a97d0b7307d4009114682f7effe2cda0f3338024bf01cae4655b7317ed4
basis_order = 2**633 * 3**6 * 5**4 * 7**3 * 11**2 * 13**2 * 17**2 * 19**2 * 29**2 * 37**2 * 41**2 * 43**2 * 59 * 71**3 * 83**2 * 89**2 * 127**2 * 139**2 * 167**2 * 211**2 * 251**2 * 281**2 * 419**2 * \
    421**2 * 503**2 * 631**2 * 701**2 * 839**2 * 1009**2 * 1259**2 * 2309**2 * 2521**2 * 2647**2 * \
    2729**2 * 3023 * 3061**2 * 3359 * 3779**2 * 4409 * \
    5039 * 6299 * 6719 * 9181 * 19531 * 22679 * 41161

# Convert a little-endian array to the number


def Array2Num(array, nbyte):
    num = 0
    for element in array[::-1]:
        num += element
        num <<= (nbyte << 3)
    num >>= (nbyte << 3)
    return num


# Convert a number to the little-endian array
def Num2Array(num, nbyte, N=0):
    array = []
    while num > 0:
        element = num & ((1 << (nbyte << 3)) - 1)
        array.append(element)
        num >>= (nbyte << 3)
    if N != 0:
        for _ in range(N - len(array)):
            array.append(0)
    return array


# Print an array in rust hex array format
def PrintArray(array, is_long=False):
    result = "["
    for element in array:
        if is_long:
            result += "0x%016X, " % (element)
        else:
            result += "0x%08X, " % (element)
    result = result[:-2] + "]"
    return result


prime = Array2Num(p_array, 8)
print("P :", hex(prime))
print("d1 :", hex(d1))
print("d2 :", hex(d2))
print("dA1 :", hex(dA1))
print("dA2 :", hex(dA2))

# print("const Px : Fq = Fq::new(%s);" %
#       (PrintArray(Num2Array(Px, 8, 21), True)))
# print("const Py : Fq = Fq::new(%s);" %
#       (PrintArray(Num2Array(Py, 8, 21), True)))
# print("const Qx : Fq = Fq::new(%s);" %
#       (PrintArray(Num2Array(Qx, 8, 21), True)))
# print("const Qy : Fq = Fq::new(%s);" %
#       (PrintArray(Num2Array(Qy, 8, 21), True)))

print(PrintArray(Num2Array(basis_order, 4)))

# PrintArray(Num2Array(d1, 4))
# PrintArray(Num2Array(d2, 4))
# PrintArray(Num2Array(dA, 4))
# PrintArray(Num2Array(dA1, 4))
# PrintArray(Num2Array(dA2, 4))
