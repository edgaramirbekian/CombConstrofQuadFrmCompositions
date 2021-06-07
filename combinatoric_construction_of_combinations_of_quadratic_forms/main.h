#include <string>
#include <vector>

using namespace std;

class CombConstrOfQuadFormComposition
{
public:
	CombConstrOfQuadFormComposition(int n_num, int r_num);
	vector<string>	CombinatorConstruction();
	static int		CalcRadonHurwitzNumber(int n);
	static bool		IsNumeric(string str);

private:
	int		n;
	int		r;
	int			SignFunction(int k, int a);
	int			TranspozitionFunction(int k, int a);
	int			CalcAuxiliaryNum_Am(int a, unsigned int m);
	int			CalculateAuxiliarySign_Tm(int a, int m);
};
