function millipos = pos_from_nibble(n0,n1,n2,n3,ar)

n0 = double(n0);
n1 = double(n1);
n2 = double(n2);
n3 = double(n3);

chomp = bitand(n0, 15) + bitshift(bitand(n1, 15), 4) + bitshift(bitand(n2,15),8) + bitshift(bitand(n3,15),12);
millipos = chomp * ar.range / ar.domain;