# This program is designed to compute the solutions for the polynomial
# equations of the degree two or less

import sys
import os
import re
import argparse

MAX_POLYNOMIAL_DEGREE = 200


def reader(equations):
	"""
		Current function checks format of input argument
		and tries to read if the file is passed
	"""
	assert len(equations) == 1, f"\033[1;91mERROR: An equation must be a single string"\
		" or a file_name with equations in a file\033[0m"
	if os.path.isfile(equations[0]):
		with open(equations[0], 'r') as input_file:
			read_file = input_file.read()
			lines = read_file.count('\n')
			equations_list = read_file.splitlines()
			print(f'{lines} equation(s) passed inside a file\n')
			return equations_list
	else:
		return equations


def normalize_signs(equation):
	"""
		For any equation -a*X^2-b*X-c=0
		reduce to normal form +-a*X^2+-b*X+-c=0
	"""
	# Regex explanation: Negative lookbehind '?<!' to match
	# '-' unless it is preceded by '^' or it's the first character
	equation = re.sub(r'(?<!\^)(?<!^)-', '+-', equation)
	return equation


def move_left_everything(equation):
	"""
		For equation with non-zero left side -a*X^2-b*X-c=-c*X^2-d*X-e
		Transform it to -a*X^2-b*X-c+c*X^2+d*X+e=0
	"""
	right_flipped = ''
	two_sides = equation.split('=')
	assert len(two_sides) == 2, """\033[1;91mERROR: Equation must have single '=' sign\033[0m"""
	assert two_sides[0] and two_sides[1], """\033[1;91mERROR: Equation left and right sides must be non-empty\033[0m"""
	if two_sides[1] != '0':
		right_flipped = two_sides[1].replace('-', '%temp%').replace('+', '-').replace('%temp%', '+')
	if right_flipped and not right_flipped[0] in ['-', '+']:
		right_flipped = '-' + right_flipped
	return two_sides[0] + right_flipped


def count_terms(equat_orig):
	"""
		Count number of terms in the original polynom
	"""
	terms = 0
	for pow_group in equat_orig:
		if pow_group:
			for _ in pow_group:
				terms += 1
	print(f'\033[1;95mTerms in the polynom: \033[0m{terms}')


def x_pow_one(nom):
	"""
		Normalize polynom by adding power value '1'
		if missing
	"""
	if 'X' in nom and '^' not in nom:
		return nom.replace('X', 'X^1')
	return nom


def x_add_one(nom):
	"""
		Normalize polynom by adding scalar '1' to 'X'
		if missing
	"""
	if 'X' in nom and '*' not in nom:
		return nom.replace('X', '1*X')
	return nom


def x_pow_zero(nom):
	"""
		Normalize polynom by adding 'X^0' to scalars
		if missing
	"""
	if 'X' not in nom:
		return nom + '*X^0'
	return nom


def get_polyn_degree(equat_simp):
	"""
		Get polynom order (degree)
	"""
	polyn_deg = 0
	for power, arg in enumerate(equat_simp):
		if arg and power > polyn_deg:
			polyn_deg = power
	return polyn_deg


def is_integer(number: float):
	"""
		Check if float can be converted to an integer
		and return either float or int
	"""
	if number.is_integer():
		return int(number)
	return number


def compute_roots(equat_simp, polyn_deg, verbose, round_to):
	"""
		Solve formated polynom
	"""
	print(f'\033[1;92mPolynomial degree: \033[0m{polyn_deg}')
	if polyn_deg > 2:
		print("""\033[1;94mThe polynomial degree is strictly greater than 2, I can't solve.\033[0m""")
	elif polyn_deg == 0:
		if equat_simp[0] == 0:
			print("""\033[1;94mAll real numbers are solutions\033[0m""")
		else:
			print(f'\033[1;94mThere is no solution\033[0m')
	elif polyn_deg == 1:
		print(f'\033[1;94mThe solution is:\033[0m\n{round(is_integer(-equat_simp[0] / equat_simp[1]), round_to)}')
	else:
		discriminant = equat_simp[1] ** 2 - 4 * equat_simp[0] * equat_simp[2]
		if verbose:
			print(f'\033[95;1mDiscriminant: \033[0m{discriminant}')
		if discriminant > 0:
			print(f'\033[1;94mDiscriminant is strictly positive, the two solutions are:\033[0m')
			x_1 = is_integer((-equat_simp[1] - discriminant**(1/2)) / (2 * equat_simp[2]))
			x_2 = is_integer((-equat_simp[1] + discriminant**(1/2)) / (2 * equat_simp[2]))
			print(round(x_1, round_to))
			print(round(x_2, round_to))
		elif discriminant == 0:
			print(f'\033[1;94mDiscriminant is zero, the single solution is:\033[0m')
			x_0 = is_integer(-equat_simp[1] / (2 * equat_simp[2]))
			print(round(x_0, round_to))
		else:
			print(f'\033[1;94mDiscriminant is strictly negative, the two solutions are complex numbers:\033[0m')
			real = is_integer(-equat_simp[1] / (2 * equat_simp[2]))
			imag = is_integer((-discriminant)**(1/2)/(2 * equat_simp[2]))
			print(f'{round(real, round_to)} - {round(imag, round_to)}i')
			print(f'{round(real, round_to)} + {round(imag, round_to)}i')

def format_equation_output(equat_simp, polyn_deg, short):
	"""
		Format equation for display in the terminal
		:short: -- <bool> to specify if terms with
		zero scalar can be omitted
	"""
	equat_format_out = ''
	for idx, coef in enumerate(equat_simp):
		if isinstance(coef, float) and coef.is_integer():
			coef = int(coef)
		if idx == 0:
			if coef < 0:
				equat_format_out = '- '
			if not coef and short:
				pass
			else:
				equat_format_out = equat_format_out + str(abs(coef)) + ' * X^0'
		elif idx > 0 and idx <= polyn_deg:
			if not coef and short:
				pass
			else:
				if coef < 0:
					if equat_format_out:
						equat_format_out = equat_format_out + ' - '
					else:
						equat_format_out = equat_format_out + '- '
				else:
					if equat_format_out:
						equat_format_out = equat_format_out + ' + '
				equat_format_out = equat_format_out + str(abs(coef)) + ' * X^' + str(idx)
	if not equat_format_out:
		equat_format_out = '0'
	return equat_format_out + ' = 0'

def normalize(equation):
	"""
		Normalize the equation:
		- remove spaces
		- move all terms to the left
		- replace sign '-' by '+-' to simplify parsing
		- add power 1 to X if it's missing
		- add X^0 to a scalar if it's missing
		- add scalar ont to X^p if it's missing
	"""
	equation = equation.replace(' ', '').upper()
	equation = move_left_everything(equation)
	equation = normalize_signs(equation)
	noms = equation.split('+')
	noms = [x_pow_one(nom) for nom in noms]
	noms = [x_pow_zero(nom) for nom in noms]
	noms = [x_add_one(nom) for nom in noms]
	return noms
   
def convert_to_list(noms):
	"""
		Convert string polynom to a list with the structure:
		list = [[term_1, term2..], [term_1, term2..], ..., [term_1, term2..]]
		where 'idx' for list[idx] corresponds to term's degree
	"""
	equat_orig = [[] for i in range(MAX_POLYNOMIAL_DEGREE + 1)]
	for nom in noms:
		nom_split = nom.split('*X^')
		assert len(nom_split) == 2, f'\033[1;91mERROR: nom {nom} must be in a format a * X^p\033[0m'
		assert len(nom_split) != 0, f'\033[1;91mERROR: there is an empty nom\033[0m'
		try:
			power = float(nom_split[1])
			if not power.is_integer():
				print(f"""\033[1;91mWrong power value <{power}> must be an integer, skipping...\033[0m""")
				return
			power = int(power)
			if power > MAX_POLYNOMIAL_DEGREE:
				print(f"\033[1;91mThe polynomial degree <{power}> must be "\
					"below {MAX_POLYNOMIAL_DEGREE}, skipping...\033[0m")
				return
			elif int(power) < 0:
				print(f"""\033[1;91mWrong power value <{power}> must be non-negative, skipping...\033[0m""")
				return
		except ValueError:
			print(f"\033[1;91mWrong power value <{nom_split[1]}> must be a number, exiting...\033[0m")
		equat_orig[int(nom_split[1])].append(nom_split[0])
	return equat_orig

def equat_reduce(equat_orig):
	"""
		Reduce polynom to the form:
		a * X^2 + b * X + c = 0
	"""
	equat_reduced = [0 for i in range(MAX_POLYNOMIAL_DEGREE + 1)]
	for idx, nom_args in enumerate(equat_orig):
		if nom_args:
			try:
				equat_reduced[idx] = sum([float(arg) for arg in nom_args])
			except ValueError:
				print(f"""\033[1;91mWrong format, equation shall contain one sign per term, exiting...\033[0m""")
				sys.exit(0)
	return equat_reduced

def solver(equation, verbose, short, round_to):
	"""
		Solver expect a polynom of noms in a format 'a * X^p'
		Where p is within [0; 2]. Spaces are arbitrary. Order is required
	"""
	noms = normalize(equation)
	equat_orig = convert_to_list(noms)
	if not equat_orig:
		return
	if verbose:
		count_terms(equat_orig)
	equat_reduced = equat_reduce(equat_orig)
	polyn_deg = get_polyn_degree(equat_reduced)
	print(f'\033[93mReduced form is: \033[0m{format_equation_output(equat_reduced, polyn_deg, short)}')
	compute_roots(equat_reduced, polyn_deg, verbose, round_to)
	return

def usage_msd():
	return ''' python3 computor.py [-v] [-s] [-r <int>] [equation|file_name]
	shere an equation or a file_name must be the last argument
	if file_name is provided -- file must contain an equation per line
	with no empty lines between equations and it must end with an empty line'''

def arg_parser():
	"""
		Parse input arguments before the equation
	"""
	parser = argparse.ArgumentParser(add_help=False, usage=usage_msd())
	parser.add_argument('-h', '--help', action='help', default=argparse.SUPPRESS,
						help='show this help')
	parser.add_argument('-v', '--verbose', help='print intermediate steps',
						action='store_true', default=False)
	parser.add_argument('-s', '--short', help='print short redused form',
						action='store_true', default=False)
	parser.add_argument('-r', '--round', help='specify N decimals to be rounded to',
						type=int, default=6)
	args, equations = parser.parse_known_args()
	if not equations:
		parser.print_help()
		sys.exit(0)
	return args, equations

def main():
	args, equations = arg_parser()
	equations = reader(equations)
	if type(equations) == list:
		for idx, equation in enumerate(equations):
			solver(equation, args.verbose, args.short, args.round)
			if idx != len(equations) - 1:
				print('')
	else:
		solver(equations, args.verbose, args.short, args.round)

if __name__ == '__main__':
	main()
