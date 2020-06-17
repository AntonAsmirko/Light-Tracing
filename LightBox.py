import math
import copy
import sys

def read_input():
    global abcd, light_vector, in_point, light_energy, num_mirrors, mirrors
    f = open('input.txt', 'r')
    inp_file = f.readlines()
    f.close()
    for i, j in zip(['a', 'b', 'c', 'd'], range(4)):
        abcd[i] = [float(i) for i in inp_file[j].split(' ')]
    light_vector = [float(i) for i in inp_file[4].split()]
    in_point = [float(i) for i in inp_file[5].split()]
    light_energy = int(inp_file[6].split()[0])
    num_mirrors = int(inp_file[7].split()[0])
    count = 9
    for i in range(num_mirrors):
        tmp = []
        for j in range(count, count + 3):
            tmp.append([float(i) for i in inp_file[j-1].split(' ')])
        mirrors.append(tmp)
        count+=3
    pass

def write_output(light_power, point):
    global light_vector
    output = open('output.txt', 'w')
    if light_power == 0:
        output.write('0\n')
        output.write(str(point[0])+' '+str(point[1])+' '+str(point[2])+'\n')
    else:
        output.write('1\n')
        output.write(str(light_power)+'\n')
        output.write(str(point[0])+' '+str(point[1])+' '+str(point[2])+'\n')
        light_vector = get_onesized(light_vector)
        output.write(str(light_vector[0])+' '+str(light_vector[1])+' '+str(light_vector[2])+'\n')
    output.close()
    pass

def build_cube_sides():
    global abcd
    floor_mirror = compute_mirror_coefficients([abcd['a'], abcd['b'], abcd['c']])
    first_mirror = compute_mirror_coefficients([abcd['b'], abcd['c'], abcd['d']])
    c_b_vector = vector_sum(abcd['b'],vector_constant_mult(abcd['c'], -1.))
    second_mirror = compute_mirror_coefficients([vector_sum(abcd['d'], c_b_vector), abcd['a'], abcd['b']])
    b_a_vector = vector_sum(abcd['a'],vector_constant_mult(abcd['b'], -1.))
    third_mirror = compute_mirror_coefficients([vector_sum(abcd['d'], b_a_vector), abcd['c'], abcd['d']])
    c_d_vector = vector_sum(abcd['d'],vector_constant_mult(abcd['c'], -1.))
    fourth_mirror = compute_mirror_coefficients([vector_sum(abcd['a'], c_d_vector), vector_sum(abcd['d'], b_a_vector), abcd['a']])
    supreme_mirror = compute_mirror_coefficients([vector_sum(abcd['a'], c_d_vector), vector_sum(abcd['d'], b_a_vector), abcd['d']])
    return [floor_mirror, first_mirror, second_mirror, third_mirror, fourth_mirror, supreme_mirror]

def vector_sum(a1, a2):
    return [a1[0]+ a2[0], a1[1]+ a2[1], a1[2]+ a2[2]]

def check_intersection(vector, plane):
    return scalar_mult(vector, [plane[0], plane[1], plane[2]]) != 0

def isclose(a, b, rel_tol=1e-04, abs_tol=0.0):
    return abs(a-b) <= max(rel_tol * max(abs(a), abs(b)), abs_tol)

def check_vector_equals(a1, a2):
    for i, j in zip(a1, a2):
        if not isclose(i, j):
            return False
    return True

def check_point_place(cur_point, start_point):
    dir_vector = vector_sum(cur_point, vector_constant_mult(start_point, -1.))
    dir_onesized = get_onesized(dir_vector)
    light_onesized = get_onesized(light_vector)
    return check_vector_equals(dir_onesized, light_onesized)

def compute_t(mirror_c):
    return -(mirror_c[0]*in_point[0]+mirror_c[1]*in_point[1]+mirror_c[2]*in_point[2]\
             +mirror_c[3])/(mirror_c[0]*light_vector[0]+mirror_c[1]*light_vector[1]+mirror_c[2]*light_vector[2])

def compute_cross_point(mirror_c):
    t = compute_t(mirror_c)
    return (in_point[0]+light_vector[0]*t, in_point[1]+light_vector[1]*t, in_point[2]+light_vector[2]*t)

def compute_mirror_coefficients(mirror):
    A = (mirror[1][1] - mirror[0][1])*(mirror[2][2] - mirror[0][2]) - (mirror[2][1] - mirror[0][1])*(mirror[1][2] - mirror[0][2])
    B = -((mirror[1][0] - mirror[0][0])*(mirror[2][2] - mirror[0][2]) - (mirror[2][0] - mirror[0][0])*(mirror[1][2] - mirror[0][2]))
    C = (mirror[1][0] - mirror[0][0])*(mirror[2][1] - mirror[0][1]) - (mirror[2][0] - mirror[0][0])*(mirror[1][1] - mirror[0][1])
    D = -mirror[0][0]*A - mirror[0][1]*B - mirror[0][2]*C
    return (A, B, C, D)



def dist(p1, p2):
    return math.sqrt((p1[0]-p2[0])**2+(p1[1]-p2[1])**2+(p1[2]-p2[2])**2)

def is_point_in_plane(plane, point):
    return plane[0]*point[0] + plane[1]*point[1] + plane[2]*point[2] + plane[3] == 0

def compute_cos(a1, a2):
    return (scalar_mult(a1, a2))/(vector_len(a1)*vector_len(a2))

def scalar_mult(a1, a2):
    return a1[0]*a2[0]+a1[1]*a2[1]+a1[2]*a2[2]

def vector_len(a1):
    return math.sqrt(a1[0]**2+a1[1]**2+a1[2]**2)

def vector_constant_mult(vector, const):
    return [i*const for i in copy.deepcopy(vector)]

def get_onesized(vector):
    return vector_constant_mult(vector, 1./vector_len(vector))

def make_projection(frm, to):
    return vector_constant_mult(to, float(scalar_mult(frm, to))/scalar_mult(to, to))

def another_reflection(plane_norm, vector):
    nrm = get_onesized(plane_norm)
    return vector_sum(vector, vector_constant_mult(nrm, 2.*scalar_mult(vector, nrm)))

def compute_reflection(plane_norm, vector, point, init_point):
    projection = vector_constant_mult(make_projection(vector_sum(point, vector_constant_mult(init_point, -1.)), plane_norm), -1.)
    proj_len = vector_len(projection)
    #onesized_norm = get_onesized(plane_norm)
    #comp_sized_norm = vector_constant_mult(onesized_norm, proj_len)
    comp_sized_norm_end_point = vector_sum(point, projection)
    proj_vector = vector_sum(comp_sized_norm_end_point, vector_constant_mult(init_point, -1.))
    proj_vector_2x = vector_constant_mult(proj_vector, 2.)
    result_end_point = vector_sum(proj_vector_2x, init_point)
    # TODO fix bag
    return vector_sum(result_end_point, vector_constant_mult(point, -1.))

def next_mirror():
    nearest_mirror_c = None
    nearest_point = None
    is_mirror = True
    min_dist = sys.maxsize
    for mirror in mirrors:
        if check_intersection(light_vector ,compute_mirror_coefficients(mirror)):
            point = compute_cross_point(compute_mirror_coefficients(mirror))
        else:
            continue
        if not is_point_in_plane(compute_mirror_coefficients(mirror), in_point) and check_point_place(point, in_point):
            dist_ = dist(point, in_point)
            if dist_ < min_dist:
                min_dist = dist_
                nearest_mirror_c = compute_mirror_coefficients(mirror)
                nearest_point = point
    cube_sides = build_cube_sides()
    for mirror in cube_sides:
        if check_intersection(light_vector, mirror):
            point = compute_cross_point(mirror)
        else:
            continue
        if not is_point_in_plane(mirror, in_point) and check_point_place(point, in_point):
            dist_ = dist(point, in_point)
            if dist_ < min_dist:
                is_mirror = False
                min_dist = dist_
                nearest_mirror_c = mirror
                nearest_point = point
    return nearest_mirror_c, nearest_point, is_mirror


abcd = {}
light_vector = [2., 4., 3.]
in_point = [5., 10., 1.]
light_energy = None
num_mirrors = None
mirrors = []


def main():
    global light_vector, light_energy, in_point, num_mirrors
    read_input()
    light_is_in_cube = True
    intersection_point = None
    while light_is_in_cube and light_energy > 0:
        nearest_mirror_c, intersection_point, light_is_in_cube = next_mirror()
        if not light_is_in_cube:
            break
        light_vector = compute_reflection([nearest_mirror_c[0], nearest_mirror_c[1], nearest_mirror_c[2]], light_vector,
                                          intersection_point, in_point)
        light_energy -= 1
        in_point = intersection_point
    write_output(light_energy, intersection_point)
    pass

if __name__ == "__main__":
    main()

#print(compute_reflection([1.,-1.,0.], [1., 0., 0.], [2.5, 0.5, 0.5],[0., 0.5, 0.5]))
#print(another_reflection(vector_constant_mult[1.,-1.,0.],-1.), [1., 0., 0.]))

#print(make_projection([0., 1., 3.], [0., 0., 1.]))