enc_str_2 = 'xKZl_^_XCY^CIE'
enc_str_1 = 'KFFSE_XHKYOKXOHOFEDM^E_Y'
byte = 42

for i in xrange(100):
    l = [chr(ord(c) ^ i) for c in enc_str_2]
    my_s = ''
    for c in l:
        my_s += c
        
    print my_s