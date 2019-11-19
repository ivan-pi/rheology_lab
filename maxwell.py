# Generalized maxwell models with varying number of elements

# 1 Maxwell element
def lG1fit(omega, g1,l1):
    return (g1*l1*omega/(1.0+l1**2*omega**2));

def sG1fit(omega, g1,l1):
    return (g1*l1**2*omega**2/(1.0+l1**2*omega**2));

def lG1eval(omega, p):
    return (p[0]*p[1]*omega/(1.0+p[1]**2*omega**2));

def sG1eval(omega, p):
    return (p[0]*p[1]**2*omega**2/(1.0+p[1]**2*omega**2));

# 2 Maxwell element
def lG2fit(omega, g1,g2,l1,l2):
    return (g1*l1*omega/(1.0+l1**2*omega**2)+
            g2*l2*omega/(1.0+l2**2*omega**2));

def sG2fit(omega, g1,g2,l1,l2):
    return (g1*l1**2*omega**2/(1.0+l1**2*omega**2)+
            g2*l2**2*omega**2/(1.0+l2**2*omega**2));

def lG2eval(omega, p):
    return (p[0]*p[2]*omega/(1.0+p[2]**2*omega**2)+
            p[1]*p[3]*omega/(1.0+p[3]**2*omega**2));

def sG2eval(omega, p):
    return (p[0]*p[2]**2*omega**2/(1.0+p[2]**2*omega**2)+
            p[1]*p[3]**2*omega**2/(1.0+p[3]**2*omega**2));

# 3 Maxwell element
def lG3fit(omega, g1,g2,g3,l1,l2,l3):
    return (g1*l1*omega/(1.0+l1**2*omega**2)+
            g2*l2*omega/(1.0+l2**2*omega**2)+
            g3*l2*omega/(1.0+l3**2*omega**2));

def sG3fit(omega, g1,g2,g3,l1,l2,l3):
    return (g1*l1**2*omega**2/(1.0+l1**2*omega**2)+
            g2*l2**2*omega**2/(1.0+l2**2*omega**2)+
            g3*l3**2*omega**2/(1.0+l3**2*omega**2));

def lG3eval(omega, p):
    return (p[0]*p[3]*omega/(1.0+p[3]**2*omega**2)+
            p[1]*p[4]*omega/(1.0+p[4]**2*omega**2)+
            p[2]*p[5]*omega/(1.0+p[5]**2*omega**2));

def sG3eval(omega, p):
    return (p[0]*p[3]**2*omega**2/(1.0+p[3]**2*omega**2)+
            p[1]*p[4]**2*omega**2/(1.0+p[4]**2*omega**2)+
            p[2]*p[5]**2*omega**2/(1.0+p[5]**2*omega**2));

# 4 Maxwell elements

def lG4fit(omega, g1,g2,g3,g4,l1,l2,l3,l4):
    return (g1*l1*omega/(1.0+l1**2*omega**2)+
            g2*l2*omega/(1.0+l2**2*omega**2)+
            g3*l3*omega/(1.0+l3**2*omega**2)+
            g4*l4*omega/(1.0+l4**2*omega**2));

def sG4fit(omega, g1,g2,g3,g4,l1,l2,l3,l4):
    return (g1*l1**2*omega**2/(1.0+l1**2*omega**2)+
            g2*l2**2*omega**2/(1.0+l2**2*omega**2)+
            g3*l3**2*omega**2/(1.0+l3**2*omega**2)+
            g4*l4**2*omega**2/(1.0+l4**2*omega**2));

def lG4eval(omega, p):
    return (p[0]*p[4]*omega/(1.0+p[4]**2*omega**2)+
            p[1]*p[5]*omega/(1.0+p[5]**2*omega**2)+
            p[2]*p[6]*omega/(1.0+p[6]**2*omega**2)+
            p[3]*p[7]*omega/(1.0+p[7]**2*omega**2));

def sG4eval(omega, p):
    return (p[0]*p[4]**2*omega**2/(1.0+p[4]**2*omega**2)+
            p[1]*p[5]**2*omega**2/(1.0+p[5]**2*omega**2)+
            p[2]*p[6]**2*omega**2/(1.0+p[6]**2*omega**2)+
            p[3]*p[7]**2*omega**2/(1.0+p[7]**2*omega**2));

# 5 Maxwell elements

def lG5fit(omega, g1,g2,g3,g4,g5,l1,l2,l3,l4,l5):
    return (g1*l1*omega/(1.0+l1**2*omega**2)+
            g2*l2*omega/(1.0+l2**2*omega**2)+
            g3*l3*omega/(1.0+l3**2*omega**2)+
            g4*l4*omega/(1.0+l4**2*omega**2)+
            g5*l5*omega/(1.0+l5**2*omega**2));

def sG5fit(omega, g1,g2,g3,g4,g5,l1,l2,l3,l4,l5):
    return (g1*l1**2*omega**2/(1.0+l1**2*omega**2)+
            g2*l2**2*omega**2/(1.0+l2**2*omega**2)+
            g3*l3**2*omega**2/(1.0+l3**2*omega**2)+
            g4*l4**2*omega**2/(1.0+l4**2*omega**2)+
            g5*l5**2*omega**2/(1.0+l5**2*omega**2));

def lG5eval(omega, p):
    return (p[0]*p[5]*omega/(1.0+p[5]**2*omega**2)+
            p[1]*p[6]*omega/(1.0+p[6]**2*omega**2)+
            p[2]*p[7]*omega/(1.0+p[7]**2*omega**2)+
            p[3]*p[8]*omega/(1.0+p[8]**2*omega**2)+
            p[4]*p[9]*omega/(1.0+p[9]**2*omega**2));

def sG5eval(omega, p):
    return (p[0]*p[5]**2*omega**2/(1.0+p[5]**2*omega**2)+
            p[1]*p[6]**2*omega**2/(1.0+p[6]**2*omega**2)+
            p[2]*p[7]**2*omega**2/(1.0+p[7]**2*omega**2)+
            p[3]*p[8]**2*omega**2/(1.0+p[8]**2*omega**2)+
            p[4]*p[9]**2*omega**2/(1.0+p[9]**2*omega**2));

# 6 Maxwell elements

def lG6fit(omega, g1,g2,g3,g4,g5,g6,l1,l2,l3,l4,l5,l6):
    return (g1*l1*omega/(1.0+l1**2*omega**2)+
            g2*l2*omega/(1.0+l2**2*omega**2)+
            g3*l3*omega/(1.0+l3**2*omega**2)+
            g4*l4*omega/(1.0+l4**2*omega**2)+
            g5*l5*omega/(1.0+l5**2*omega**2)+
            g6*l6*omega/(1.0+l6**2*omega**2));

def sG6fit(omega, g1,g2,g3,g4,g5,g6,l1,l2,l3,l4,l5,l6):
    return (g1*l1**2*omega**2/(1.0+l1**2*omega**2)+
            g2*l2**2*omega**2/(1.0+l2**2*omega**2)+
            g3*l3**2*omega**2/(1.0+l3**2*omega**2)+
            g4*l4**2*omega**2/(1.0+l4**2*omega**2)+
            g5*l5**2*omega**2/(1.0+l5**2*omega**2)+
            g6*l6**2*omega**2/(1.0+l6**2*omega**2));

def lG6eval(omega, p):
    return (p[0]*p[6]*omega/(1.0+p[6]**2*omega**2)+
            p[1]*p[7]*omega/(1.0+p[7]**2*omega**2)+
            p[2]*p[8]*omega/(1.0+p[8]**2*omega**2)+
            p[3]*p[9]*omega/(1.0+p[9]**2*omega**2)+
            p[4]*p[10]*omega/(1.0+p[10]**2*omega**2)+
            p[5]*p[11]*omega/(1.0+p[11]**2*omega**2));

def sG6eval(omega, p):
    return (p[0]*p[6]**2*omega**2/(1.0+p[6]**2*omega**2)+
            p[1]*p[7]**2*omega**2/(1.0+p[7]**2*omega**2)+
            p[2]*p[8]**2*omega**2/(1.0+p[8]**2*omega**2)+
            p[3]*p[9]**2*omega**2/(1.0+p[9]**2*omega**2)+
            p[4]*p[10]**2*omega**2/(1.0+p[10]**2*omega**2)+
            p[5]*p[11]**2*omega**2/(1.0+p[11]**2*omega**2));

# 7 Maxwell elements

def lG7fit(omega, g1,g2,g3,g4,g5,g6,g7,l1,l2,l3,l4,l5,l6,l7):
    return (g1*l1*omega/(1.0+l1**2*omega**2)+
            g2*l2*omega/(1.0+l2**2*omega**2)+
            g3*l3*omega/(1.0+l3**2*omega**2)+
            g4*l4*omega/(1.0+l4**2*omega**2)+
            g5*l5*omega/(1.0+l5**2*omega**2)+
            g6*l6*omega/(1.0+l6**2*omega**2)+
            g7*l7*omega/(1.0+l7**2*omega**2));

def sG7fit(omega, g1,g2,g3,g4,g5,g6,g7,l1,l2,l3,l4,l5,l6,l7):
    return (g1*l1**2*omega**2/(1.0+l1**2*omega**2)+
            g2*l2**2*omega**2/(1.0+l2**2*omega**2)+
            g3*l3**2*omega**2/(1.0+l3**2*omega**2)+
            g4*l4**2*omega**2/(1.0+l4**2*omega**2)+
            g5*l5**2*omega**2/(1.0+l5**2*omega**2)+
            g6*l6**2*omega**2/(1.0+l6**2*omega**2)+
            g7*l7**2*omega**2/(1.0+l7**2*omega**2));

def lG7eval(omega, p):
    return (p[0]*p[7]*omega/(1.0+p[7]**2*omega**2)+
            p[1]*p[8]*omega/(1.0+p[8]**2*omega**2)+
            p[2]*p[9]*omega/(1.0+p[9]**2*omega**2)+
            p[3]*p[10]*omega/(1.0+p[10]**2*omega**2)+
            p[4]*p[11]*omega/(1.0+p[11]**2*omega**2)+
            p[5]*p[12]*omega/(1.0+p[12]**2*omega**2)+
            p[6]*p[13]*omega/(1.0+p[13]**2*omega**2));

def sG7eval(omega, p):
    return (p[0]*p[7]**2*omega**2/(1.0+p[7]**2*omega**2)+
            p[1]*p[8]**2*omega**2/(1.0+p[8]**2*omega**2)+
            p[2]*p[9]**2*omega**2/(1.0+p[9]**2*omega**2)+
            p[3]*p[10]**2*omega**2/(1.0+p[10]**2*omega**2)+
            p[4]*p[11]**2*omega**2/(1.0+p[11]**2*omega**2)+
            p[5]*p[12]**2*omega**2/(1.0+p[12]**2*omega**2)+
            p[6]*p[13]**2*omega**2/(1.0+p[13]**2*omega**2));

# 8 Maxwell elements

def lG8fit(omega, g1,g2,g3,g4,g5,g6,g7,g8,l1,l2,l3,l4,l5,l6,l7,l8):
    return (g1*l1*omega/(1.0+l1**2*omega**2)+
            g2*l2*omega/(1.0+l2**2*omega**2)+
            g3*l3*omega/(1.0+l3**2*omega**2)+
            g4*l4*omega/(1.0+l4**2*omega**2)+
            g5*l5*omega/(1.0+l5**2*omega**2)+
            g6*l6*omega/(1.0+l6**2*omega**2)+
            g7*l7*omega/(1.0+l7**2*omega**2)+
            g8*l8*omega/(1.0+l8**2*omega**2));

def sG8fit(omega, g1,g2,g3,g4,g5,g6,g7,g8,l1,l2,l3,l4,l5,l6,l7,l8):
    return (g1*l1**2*omega**2/(1.0+l1**2*omega**2)+
            g2*l2**2*omega**2/(1.0+l2**2*omega**2)+
            g3*l3**2*omega**2/(1.0+l3**2*omega**2)+
            g4*l4**2*omega**2/(1.0+l4**2*omega**2)+
            g5*l5**2*omega**2/(1.0+l5**2*omega**2)+
            g6*l6**2*omega**2/(1.0+l6**2*omega**2)+
            g7*l7**2*omega**2/(1.0+l7**2*omega**2)+
            g8*l8**2*omega**2/(1.0+l8**2*omega**2));

def lG8eval(omega, p):
    return (p[0]*p[8] *omega**2/(1.0+p[8]**2*omega**2)+
            p[1]*p[9] *omega**2/(1.0+p[9]**2*omega**2)+
            p[2]*p[10]*omega**2/(1.0+p[10]**2*omega**2)+
            p[3]*p[11]*omega**2/(1.0+p[11]**2*omega**2)+
            p[4]*p[12]*omega**2/(1.0+p[12]**2*omega**2)+
            p[5]*p[13]*omega**2/(1.0+p[13]**2*omega**2)+
            p[6]*p[14]*omega**2/(1.0+p[14]**2*omega**2)+
            p[7]*p[15]*omega**2/(1.0+p[15]**2*omega**2));

def sG8eval(omega, p):
    return (p[0]*p[8]**2 *omega**2/(1.0+p[8]**2 *omega**2)+
            p[1]*p[9]**2 *omega**2/(1.0+p[9]**2 *omega**2)+
            p[2]*p[10]**2*omega**2/(1.0+p[10]**2*omega**2)+
            p[3]*p[11]**2*omega**2/(1.0+p[11]**2*omega**2)+
            p[4]*p[12]**2*omega**2/(1.0+p[12]**2*omega**2)+
            p[5]*p[13]**2*omega**2/(1.0+p[13]**2*omega**2)+
            p[6]*p[14]**2*omega**2/(1.0+p[14]**2*omega**2)+
            p[7]*p[15]**2*omega**2/(1.0+p[15]**2*omega**2));

