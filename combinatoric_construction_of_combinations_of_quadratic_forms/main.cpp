#include <bitset>
#include <cmath>
#include <iostream>
#include "main.h"

CombConstrOfQuadFormComposition::CombConstrOfQuadFormComposition(int n_num, int r_num)
{
	n = n_num - 1;
	r = r_num - 1;
}

vector<string> CombConstrOfQuadFormComposition::CombinatorConstruction()
{
	vector<string> quadCompositions;
	// i-ի 0-ից n արժեքների համար
	for (int i = 0; i <= n; ++i)
	{
		string tmpQuadComposition = "";
		// k-ի 0-ից r արժեքների համար
		for (int k = 0; k <= r; ++k)
		{
			int s = SignFunction(k, i);
			string s_str = "";
			int p = TranspozitionFunction(k, i);
			if (s == 1)
			{
				if (k != 0)
					s_str = " + ";
				else
					s_str = "";
			}
			else if (s == -1)
				s_str = " - ";
			else
				s_str = to_string(s);

			if (k == 0)
				tmpQuadComposition += "z[" + to_string(i) + "] (x,y) = ";

			tmpQuadComposition += s_str + "x[" + to_string(p) + "]" + "y[" + to_string(k) + "]";
		}
		quadCompositions.push_back(tmpQuadComposition);
	}

	return quadCompositions;
}

int CombConstrOfQuadFormComposition::CalcRadonHurwitzNumber(int n)
{
	int RadonHurwitzNum = 0;

	// կենտ n թվերի համար Ռադոն-Հուրվիցի թիվը հավասար է 1-ի
	if (n % 2 == 1)
		RadonHurwitzNum = 1;
	else
	{
		// n-ը n = 2^b * (2a + 1) տեսքով նկարագրելու համար սահմանենք
		int b, a;
		// b-ն b = 4c + d && 0 <= d <= 3 տեսքով նկարագրելու համար սահմանենք
		int c, d;

		// b-ն գտնելու համար պետք է գտնենք n-ի ամենամեծ 2-ի աստիճան բաժանարարը
		// մենք չենք օգտվի հերթականությամբ տարբերակները փորձելու մեթոդը այլ կօգտվենք բիտային հնարքից
		// եթե դիտարկենք n թիվը երկուական տեսքով ապա կտեսնենք որ մեզ պետք է գտնել այն թիվը որը ստացվում է
		// n[2-ական]-ի միայն ամենաաջ 1 բիթը պահպանելով և մյուս բոլոր բիտերը 0-ացնելով (օր․ 5=101 ամենամեծ 2-ի էքսպոնենտ)
		// բաժանարարը 2^0=1=001, 48=110000 ամենամեծ 2-ի էքսպոնենտ բաժանարարը 2^4=16=010000 և այլն)։ Կկիրառենք հետևյալ ալգորիթմը։
		// գտնենք (n-1)-ը, ~ բիտային Ժխտում օպերատորով շրջենք ստացված (n-1)-ի բիթերը, այնուհետև n և (n-1)-ի միջև
		// կիրառենք & բիտային Եւ օպերացիան և կստանանք n թվի ամենամեծ 2-ի էքպոնենտ բաժանարարը և կորոշենք թե դա 2-ի որ աստիճանն է։
		int grtstDivsrOf2Exponent = n & (~(n - 1));

		b = static_cast<int>(log2(grtstDivsrOf2Exponent));
		a = ((n / grtstDivsrOf2Exponent) - 1) / 2;

		// b-ի միջոցով գտնենք c և d
		if (b < 4)
		{
			c = 0;
			d = b;
		}
		else
		{
			d = b % 4;
			c = b / 4;
		}

		// ստացված թվերով և RandonHurwitzNum(n) = 2^d + 8c բանաձևով գտնենք ՌանդոնՀուրվիցի թիվը
		RadonHurwitzNum = static_cast<int>(pow(2, d)) + (8 * c);
	}

	return RadonHurwitzNum;
}

int CombConstrOfQuadFormComposition::SignFunction(int k, int a)
{
	vector<vector<int>> signMatrix =
	{
		{1, 1, 1, 1, 1, 1, 1, 1},
		{1, -1, 1, -1, 1, -1, -1, 1},
		{1, -1, -1, 1, 1, 1, -1, -1},
		{1, 1, -1, -1, 1, -1, 1, -1},
		{1, -1, -1, -1, -1, 1, 1, 1},
		{1, 1, -1, 1, -1, -1, -1, 1},
		{1, 1, 1, -1, -1, 1, -1, -1},
		{1, -1, 1, 1, -1, -1, 1, -1}
	};
	int S = 1;
	int m = k / 8;
	int l = k % 8;

	if (k == 0)
	{
		S = 1;
	}
	else if (m > 0 && l == 0)
	{
		S = CalculateAuxiliarySign_Tm(a, m - 1);
	}
	else if (m >= 0 && l > 0 && l < static_cast<int>(signMatrix.size()))
	{
		S = CalculateAuxiliarySign_Tm(a, m) * signMatrix[l][CalcAuxiliaryNum_Am(a, m)];
	}

	return S;
}

int CombConstrOfQuadFormComposition::TranspozitionFunction(int k, int a)
{
	string binary_k_str = bitset<32>(k).to_string();
	string binary_a_str = bitset<32>(a).to_string();
	unsigned long long int binary_k = stoull(binary_k_str);
	unsigned long long int binary_a = stoull(binary_a_str);

	unsigned long long int maxBinNum = (binary_k > binary_a) ? binary_k : binary_a;

	int P = 0;
	int it = 0;

	while (maxBinNum != 0)
	{
		// Բաժանել ստացված երկուական թիվը երկուական 2-ի (2[երկուական] = 10)
		maxBinNum /= 10;
		if (binary_k % 10 + binary_a % 10 == 1)
		{
			P += static_cast<int>(pow(2, it));
		}
		binary_k = (binary_k - (binary_k % 10)) / 10;
		binary_a = (binary_a - (binary_a % 10)) / 10;
		++it;
	}

	return P;
}

int CombConstrOfQuadFormComposition::CalcAuxiliaryNum_Am(int a, unsigned int m)
{
	// 2-ական a-ն տողից դարձնել թիվ որպեսզի ավելորդ 0-ները վերանան
	string binary_a_str = to_string(stoull(bitset<32>(a).to_string()));
	vector<int> reverse_bin_a_digits;
	for (string::const_reverse_iterator crit = binary_a_str.rbegin(); crit != binary_a_str.rend(); ++crit)
	{
		int a_digit = stoi(string(1, *crit));
		reverse_bin_a_digits.push_back(a_digit);
	}

	// գտնենք a(m)-ը հետևյալ բանաձևով a(m) = a[4m] + 2*a[4m+1] + 4*a[4m+2], m = 0, 1, ..., reverse_bin_a_digits.size()
	int am = 0;
	unsigned int idx = 4 * m;
	int idx_coefficent = 1;
	while (idx_coefficent < 5)
	{
		if (idx < reverse_bin_a_digits.size())
		{
			am += idx_coefficent * reverse_bin_a_digits[idx];
		}
		idx_coefficent *= 2;
		++idx;
	}
	return am;
}

int CombConstrOfQuadFormComposition::CalculateAuxiliarySign_Tm(int a, int m)
{
	string binary_a_str = to_string(stoull(bitset<32>(a).to_string()));
	vector<int> reverse_bin_a_digits;
	for (string::const_reverse_iterator crit = binary_a_str.rbegin(); crit != binary_a_str.rend(); ++crit)
	{
		int a_digit = stoi(string(1, *crit));
		reverse_bin_a_digits.push_back(a_digit);
	}

	int sum = 0;
	int T = 1;
	for (int i = m + 1; 4 * i - 1 < static_cast<int>(reverse_bin_a_digits.size()); i += 4)
	{
		sum += reverse_bin_a_digits[4 * i - 1];
	}

	if (sum % 2 != 0)
	{
		T = -1;
	}
	return T;
}

bool CombConstrOfQuadFormComposition::IsNumeric(string str) {
	for (unsigned int i = 0; i < str.length(); i++)
	{
		if (isdigit(str[i]) == false)
			return false;
	}
	return true;
}

int MainExecution()
{
	cout << "Enter non-negative non-zero n number from Z{1,2,...} set" << endl;
	string n_str;
	cout << "n: ";
	cin >> n_str;

	if (!CombConstrOfQuadFormComposition::IsNumeric(n_str))
	{
		cout << "The entered value for n is not numeric. Please enter only numeric characters" << endl;
		return -2;
	}

	int n = stoi(n_str);

	if (n <= 0)
	{
		cout << "The entered number n does not meet the above requirements" << endl;
		return -1;
	}
	int RadHurwNum = CombConstrOfQuadFormComposition::CalcRadonHurwitzNumber(n);
	cout << "Enter non-negative r number from Z{0, 1, ..., RadonHurwitz(n)} set " << " which is <= than " << " Radon-Hurwitz number for n (which is equal to " << RadHurwNum << ")" << endl;
	string r_str;
	cout << "r: ";
	cin >> r_str;

	if (!CombConstrOfQuadFormComposition::IsNumeric(r_str))
	{
		cout << "The entered value for r is not numeric. Please enter only numeric characters" << endl;
		return -2;
	}

	int r = stoi(r_str);

	if (r > RadHurwNum)
	{
		cout << "The entered number r does not meet the above requirements" << endl;
		return -1;
	}

	CombConstrOfQuadFormComposition combConstructor = CombConstrOfQuadFormComposition(n, r);
	vector<string> quadCompositions = combConstructor.CombinatorConstruction();
	cout << "(n = " << n << ", r = " << r << ")," << " combinatorially built z[i] : R^(n + 1) * R^(n + 1) -> R^(n + 1) type quadratic form compositions for 0 <= r <= RH(n) satisfying (n, r) number couple \n" << endl;
	for (vector<string>::const_iterator cit = quadCompositions.begin(); cit != quadCompositions.end(); ++cit)
	{
		string quadComp = *cit;
		cout << quadComp << "\n" << endl;
	}
	return 0;
}

int main(int argc, char const* argv[])
{
	int retVal = 1;
	bool continueExecution = false;
	string isUserContinue = "n";
	do
	{
		retVal = MainExecution();
		cout << "Do you want to continue? (enter y for yes, n for no): ";
		cin >> isUserContinue;
		cout << "\n";

		if (isUserContinue == "y")
			continueExecution = true;
		else
			continueExecution = false;

	} while (continueExecution);

	return 0;
}
