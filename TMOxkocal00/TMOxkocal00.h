#include "TMO.h"

class TMOxkocal00 : public TMO
{
public:
	TMOxkocal00();
	virtual ~TMOxkocal00();
	virtual int Transform();

protected:
	TMODouble dParameter;
	TMOInt dNeighbourhood;
	TMOInt dColors;
	TMOInt dBins;
	TMODouble dWa;
	TMODouble dWb;
	TMODouble dLambda;
	TMOInt dJustRecolor;
	TMOInt dNormalize;

public:
	// typedef struct rgb t_rgb;
	typedef struct rgb
	{
		double r;
		double g;
		double b;

		bool operator==(const struct rgb &o) const{
				return r == o.r && g == o.g && b == o.b;
		}

		bool operator<(const struct rgb &o) const{
			return (r + g + b) < (o.r + o.g + o.b);
		}
	} t_rgb;

	// typedef struct rgb_int t_rgb_int;
	typedef struct rgb_int
	{
		int r;
		int g;
		int b;

		bool operator==(const struct rgb_int &o) const{
				return r == o.r && g == o.g && b == o.b;
		}

		bool operator<(const struct rgb_int &o) const{
			return (r + g + b) < (o.r + o.g + o.b);
		}
	} t_rgb_int;

	typedef struct Lab t_Lab;
	// typedef struct Lab
	// {
	// 	double L;
	// 	double a;
	// 	double b;

	// 	bool operator<(const struct Lab &o) const{
	// 		return (L + a + b) < (o.L + o.a + o.b);
	// 	}
	// } t_Lab;

	typedef struct Lab t_Lab_int;
	// typedef struct Lab_int
	// {
	// 	int L;
	// 	int a;
	// 	int b;

	// 	bool operator==(const struct Lab_int &o) const{
	// 			return L == o.L && a == o.a && b == o.b;
	// 	}

	// 	bool operator<(const struct Lab_int &o) const{
	// 		return (L + a + b) < (o.L + o.a + o.b);
	// 	}
	// } t_Lab_int;

private:

	void getPixelByOffset(int y, int x, int dy, int dx, double *pixel);

	double LabColorDistance(double L1, double a1, double b1, double L2, double a2, double b2);
	double LabColorDistance(t_Lab c1, t_Lab c2);
	// void getMinMaxLab(double *img, int width, int height, double *minA, double *minB, double *maxA, double *maxB);

	int getClosestColorRGB(double r, double g, double b, t_rgb* color_space, int color_space_size);
	int getClosestColorLab(double L, double a, double b, t_Lab* color_space, int color_space_size);
	// vector<TMOxkocal00::t_rgb> populirizeColors(double *img, int height, int width, int bins, int maxColors);
	// vector<t_Lab> popularizeColorsLab(double *img, int height, int width, int bins, int maxColors);

	double deltaPrimeIJ(t_Lab i_color, t_Lab j_color, t_Lab *landmarks, double normalizer);
};
