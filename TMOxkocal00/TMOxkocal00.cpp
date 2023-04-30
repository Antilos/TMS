/* --------------------------------------------------------------------------- *
 * TMOxkocal00.cpp: implementation of the TMOYourOperatorName class.   *
 * --------------------------------------------------------------------------- */

#include "TMOxkocal00.h"
#include <cmath>
#include <vector>
#include <utility>
#include <map>
#include <algorithm>
#include <limits>
#include <boost/multi_array.hpp>
#include "laplace.h"

/* --------------------------------------------------------------------------- *
 * Constructor serves for describing a technique and input parameters          *
 * --------------------------------------------------------------------------- */
TMOxkocal00::TMOxkocal00()
{
	SetName(L"TMOxkocal00");					  // TODO - Insert operator name
	SetDescription(L"Add your TMO description here"); // TODO - Insert description

	dParameter.SetName(L"ParameterName");				// TODO - Insert parameters names
	dParameter.SetDescription(L"ParameterDescription"); // TODO - Insert parameter descriptions
	dParameter.SetDefault(1);							// TODO - Add default values
	dParameter = 1.;
	dParameter.SetRange(-1000.0, 1000.0); // TODO - Add acceptable range if needed
	this->Register(dParameter);

   dNeighbourhood.SetName(L"Neighbourhood");
   dNeighbourhood.SetDescription(L"Size of the neighbourhood");
   dNeighbourhood.SetDefault(4);
   dNeighbourhood = 4;
   dNeighbourhood.SetRange(0, 100);
   this->Register(dNeighbourhood);

   dColors.SetName(L"Colors");
   dColors.SetDescription(L"Size of the popular color space");
   dColors.SetDefault(30);
   dColors = 30;
   dColors.SetRange(0, 1000000000);
   this->Register(dColors);

   dBins.SetName(L"Bins");
   dBins.SetDescription(L"Number of bins for the popular color space");
   dBins.SetDefault(3);
   dBins = 3;
   dBins.SetRange(0, 1000000000);
   this->Register(dBins);

   dWa.SetName(L"Wa");
   dWa.SetDescription(L"Control parameter for contribution of 'a' channel to color distance");
   dWa.SetDefault(1.0);
   dWa = 1.0;
   dWa.SetRange(0, 1.0);
   this->Register(dWa);

   dWb.SetName(L"Wb");
   dWb.SetDescription(L"Control parameter for contribution of 'b' channel to color distance");
   dWb.SetDefault(1.0);
   dWb = 1.0;
   dWb.SetRange(0, 1.0);
   this->Register(dWb);

   dLambda.SetName(L"Lambda");
   dLambda.SetDescription(L"Control parameter for contribution of global contrast");
   dLambda.SetDefault(0.5);
   dLambda = 0.5;
   dLambda.SetRange(0, 0.999999999999);
   this->Register(dLambda);

   dJustRecolor.SetName(L"JustRecolor");
   dJustRecolor.SetDescription(L"Produces a new image in the reduced color space");
   dJustRecolor.SetDefault(0);
   dJustRecolor = 0;
   dJustRecolor.SetRange(0, 1);
   this->Register(dJustRecolor);

   dNormalize.SetName(L"Normalize");
   dNormalize.SetDescription(L"Whether to normalize the output");
   dNormalize.SetDefault(1);
   dNormalize = 1;
   dNormalize.SetRange(0, 2);
   this->Register(dNormalize);
}

using namespace std;

TMOxkocal00::~TMOxkocal00()
{
}

struct Lab_int
{
   int L;
   int a;
   int b;

   bool operator==(const struct Lab_int &o) const{
         return L == o.L && a == o.a && b == o.b;
   }

   bool operator<(const struct Lab_int &o) const{
      return (L + a + b) < (o.L + o.a + o.b);
   }
};

struct Lab
{
   double L;
   double a;
   double b;

   bool operator<(const struct Lab &o) const{
      return (L + a + b) < (o.L + o.a + o.b);
   }
};


double colorDistance(double r1, double g1, double b1, double r2, double g2, double b2){
   return sqrt((r2 - r1)*(r2 - r1) + (g2 - g1)*(g2 - g1) + (b2 - b1)*(b2 - b1));
}

double TMOxkocal00::LabColorDistance(double L1, double a1, double b1, double L2, double a2, double b2){
   return cbrt(pow((L2-L1), 3) + dWa * pow((a2-a1), 3) + dWb * pow((b2-b1), 3));
}

double TMOxkocal00::LabColorDistance(t_Lab c1, t_Lab c2){
   double tmp1 = pow((c2.L - c1.L), 3) + dWa * pow((c2.a - c1.a), 3) + dWb * pow((c2.b - c1.b), 3);
   double tmp2 = cbrt(tmp1);
   // cerr << "[DEBUG]" << "c1=(" <<c1.L << "," << c1.a << "," << c1.b << ")" << "c2=(" <<c2.L << "," << c2.a << "," << c2.b << ")" << "tmp1=" << tmp1 <<  " tmp2=" << tmp2 << endl; 
   return tmp2;
}

int TMOxkocal00::getClosestColorRGB(double r, double g, double b, t_rgb* color_space, int color_space_size){
   int closest_index = 0;
   double closest_distance = numeric_limits<double>::max();
   for (int i = 0; i < color_space_size; i++){
      double r2 = color_space[i].r;
      double g2 = color_space[i].g;
      double b2 = color_space[i].b;

      double new_distance = colorDistance(r, g, b, r2, g2, b2);
      if(new_distance < closest_distance){
         closest_distance = new_distance;
         closest_index = i;
      }
   }

   return closest_index;
}

int TMOxkocal00::getClosestColorLab(double L, double a, double b, t_Lab* color_space, int color_space_size){
   int closest_index = 0;
   double closest_distance = numeric_limits<double>::max();
   for (int i = 0; i < color_space_size; i++){
      double L2 = color_space[i].L;
      double a2 = color_space[i].a;
      double b2 = color_space[i].b;

      double new_distance = LabColorDistance(L, a, b, L2, a2, b2);
      if(new_distance < closest_distance){
         closest_distance = new_distance;
         closest_index = i;
      }
   }

   return closest_index;
}

void getDataByCoords(int row, int col, double *data, int width, double *first, double *second, double *third){
   *first = data[((row*3)*width + (col*3))];
   *second = data[((row*3)*width + (col*3) + 1)];
   *third = data[((row*3)*width + (col*3) + 2)];
}

void setDataByCoords(int row, int col, double *data, int width, double first, double second, double third){
   data[((row*3)*width + (col*3))] = first;
   data[((row*3)*width + (col*3) + 1)] = second;
   data[((row*3)*width + (col*3) + 2)] = third;
}

int getBin(double val, vector<double> bins){
   for (size_t i = 0; i < bins.size(); i++)
   {
      if (val < bins[i]){
         return i;
      }
   }
   return bins.size();
}

double getBinBound(int bin, vector<double> bin_bounds){
   return 0.0;
}

void getMinMaxLab(double *img, int width, int height, double *minL, double *minA, double *minB, double *maxL, double *maxA, double *maxB){
   double cur_minL, cur_minA, cur_minB = numeric_limits<double>::max();
   double cur_maxL, cur_maxA, cur_maxB = numeric_limits<double>::min();
   for (int row = 0; row < height; row++)
	{
		for (int col = 0; col < width; col++)
		{
         double L, a, b;
         getDataByCoords(row, col, img, width, &L, &a, &b);
         if (a < cur_minA)
            cur_minA = a;
         if (b < cur_minB)
            cur_minB = b;
         if (L < cur_minL)
            cur_minL = L;
         if (a > cur_maxA)
            cur_maxA = a;
         if (b > cur_maxB)
            cur_maxB = b;
         if (L > cur_maxB)
            cur_maxL = L;
      }
   }
   // *minL = cur_minL < 0.0 ? 0.0 : cur_minL;
   *minL = cur_minL;
   *minA = cur_minA;
   *minB = cur_minB;
   *maxL = cur_maxL > 100.0 ? 100.0 : cur_maxL;
   *maxA = cur_maxA;
   *maxB = cur_maxB;
}

vector<TMOxkocal00::t_Lab> popularizeColorsLab(double *img, int height, int width, int bins, int maxColors){
   //compute min and max values
   double minL, minA, minB, maxL, maxA, maxB;
   // minL = 0.0;
   // maxL = 100.0;
   getMinMaxLab(img, width, height, &minL, &minA, &minB, &maxL, &maxA, &maxB);
   cerr << minL << " " << minA << " " << minB << " " << maxL << " " << maxA << " " << maxB << " " << endl;

   //construct bins
   bins = pow(2, bins+1);
   vector<double> L_bin_bounds, A_bin_bounds, B_bin_bounds;
   double L_bin_size = (maxL - minL) / double(bins);
   double A_bin_size = (maxA - minA) / double(bins);
   double B_bin_size = (maxB - minB) / double(bins);
   // cerr << bin_size << endl;
   for (size_t i = 1; i < bins; i++)
   {
      L_bin_bounds.push_back(minL + i*L_bin_size);
      A_bin_bounds.push_back(minA + i*A_bin_size);
      B_bin_bounds.push_back(minB + i*B_bin_size);
   }
   L_bin_bounds.push_back(maxL);
   A_bin_bounds.push_back(maxA);
   B_bin_bounds.push_back(maxB);

   // cerr << "[DEBUG] Bounds: ";
   // for (size_t i = 1; i < bins; i++)
   // {
   //    cerr << bin_bounds[i] << " ";
   // }
   // cerr << endl;

   map<TMOxkocal00::t_Lab_int, int> colors_map;

   //scan image for colors
   int color_counter = 0;
   for (int j = 0; j < height; j++)
	{
		for (int i = 0; i < width; i++)
		{
         TMOxkocal00::t_Lab_int color;
         color.L = getBin(*img++, L_bin_bounds);
         color.a = getBin(*img++, A_bin_bounds);
         color.b = getBin(*img++, B_bin_bounds);

         // cerr << "[DEBUG] colorRaw=" << color.L << "(" << img[i] << ") " << color.a << "(" << img[i] << ") " << color.b << "(" << img[i] << ")" << endl;

         auto search = colors_map.find(color);
         // cerr << "DING 1" << endl;
         if(search != colors_map.end()){
            (search->second)++;
            // cerr << "DING 1.0" << endl;
         }else{
            colors_map.insert({color, 1});
            // cerr << "DING 1.5" << endl;
            color_counter++;
         }
      }
   }

   cerr << "DING 1: colors colected=" << colors_map.size() << endl;
   cerr << "DING 1: colors encountered=" << color_counter << endl;

   //sort colors
   vector<pair<TMOxkocal00::t_Lab_int, int>> colors_vector;
   for (auto& it : colors_map){
      colors_vector.push_back(it);
   }

   sort(colors_vector.begin(), colors_vector.end(), [](pair<TMOxkocal00::t_Lab_int, int>& a, pair<TMOxkocal00::t_Lab_int, int>& b){return a.second > b.second;});

   cerr << "DING 2" << endl;

   //get new color space
   vector<TMOxkocal00::t_Lab> new_colors;
   maxColors = maxColors <= colors_vector.size() ? maxColors : colors_vector.size();
   for(int i = 0; i < maxColors; i++){
      TMOxkocal00::t_Lab new_color;
      new_color.L = L_bin_bounds[colors_vector[i].first.L];
      new_color.a = A_bin_bounds[colors_vector[i].first.a];
      new_color.b = B_bin_bounds[colors_vector[i].first.b];
      new_colors.push_back(new_color);
      cerr << "[DEBUG] color=(" << new_color.L << ", " << new_color.a << ", " << new_color.b << ") pop=" << colors_vector[i].second <<endl;
   }

   cerr << "DING 3" << endl;

   return new_colors;
}

double TMOxkocal00::deltaPrimeIJ(t_Lab i_color, t_Lab j_color, t_Lab *landmarks, double normalizer){
   //calculate sum of distance from landmarks
   double landmark_distance_sum = 0.0;
   for (size_t i = 0; i < dColors; i++)
   {
      double tmp1 = LabColorDistance(i_color, landmarks[i]);
      double tmp2 = LabColorDistance(landmarks[i], j_color);
      // cerr << tmp1 << " " << tmp2 << endl;
      landmark_distance_sum += tmp1 + tmp2;
   }
   // cerr << "[DEBUG] l_distance_sum=" <<landmark_distance_sum << endl;
   

   return (((1.0 - dLambda) * LabColorDistance(i_color, j_color)) + ((dLambda / 2.0) * landmark_distance_sum)) / normalizer;
}

double getMeanL(double* img, int width, int height){
   double sum = 0.0;
   for (size_t i = 0; i < width; i++)
   {
      for (size_t j = 0; j < height; j++)
      {
         double L, a, b;
         getDataByCoords(i, j, img, width, &L, &a, &b);
         sum += L;
      }
   }
   
   return sum / (double)(width * height);
}



/* --------------------------------------------------------------------------- *
 * This overloaded function is an implementation of your tone mapping operator *
 * --------------------------------------------------------------------------- */
int TMOxkocal00::Transform()
{
   cout << "[DEBUG] Starting" << endl;

   //Calculate some constant values
   double delta_prime_normalizer = (1.0-dLambda) + ((double)(dLambda * dColors) / 2.0);
   cerr << "[DEBUG] normalizer=" << delta_prime_normalizer << endl;


	// Source image is stored in local parameter pSrc
	// Destination image is in pDst

	// Initialy images are in RGB format, but you can
	// convert it into other format
	pSrc->Convert(TMO_LAB); 
	pDst->Convert(TMO_LAB);

	double *pSourceData = pSrc->GetData();		// You can work at low level data
	double *pDestinationData = pDst->GetData(); // Data are stored in form of array
												// of three doubles representing
												// three colour components
	double pY, px, py;
   

   cerr << "[DEBUG] Starting" << endl;

   vector<t_Lab> new_colors = popularizeColorsLab(pSourceData, pSrc->GetHeight(), pSrc->GetWidth(), dBins, dColors);
   // for(t_Lab color : new_colors){
   //    cerr << "[DEBUG] color=(" << color.L << ", " << color.a << ", " << color.b << ")" <<endl;
   // }
   cerr << new_colors.size() << endl;

   cerr << "[DEBUG] Color space calculated" << endl;

   boost::multi_array<double, 2> U;
   double U_min = numeric_limits<double>::max();
   double U_max = numeric_limits<double>::min();

   // Calculate b' vector
   if (!dJustRecolor){
      vector<double> bPrime;
      boost::multi_array<double, 2> F(boost::extents[pSrc->GetHeight()][pSrc->GetWidth()]);
      cerr << "[DEBUG] starting to collect b' vector" << endl;
      for (int row = 0; row < pSrc->GetHeight(); row++)
      {
         for (int col = 0; col < pSrc->GetWidth(); col++)
         {
            double L, a, b;
            getDataByCoords(row, col, pSourceData, pSrc->GetWidth(), &L, &a, &b);
            t_Lab color;
            color.L = L;
            color.a = a;
            color.b = b;

            //colect neighbours
            double neigh_distances = 0.0;
            for (int ii = -1; ii <= 1; ii++){
               for (int jj = -1; jj <= 1; jj++){
                  if (abs(ii) == abs(jj)){
                     continue;
                  }
                  
                  if (row + ii < 0 || row + ii >= pSrc->GetHeight() || col+jj < 0 || col+jj >= pSrc->GetWidth()){
                     continue;
                  }


                  double Lj, aj, bj;
                  // cerr << "Getting neighbour" << endl;
                  getDataByCoords(row+ii, col+jj, pSourceData, pSrc->GetWidth(), &Lj, &aj, &bj);
                  // cerr << "Ha got the neighbour" << endl;
                  t_Lab neigh;
                  neigh.L = Lj;
                  neigh.a = aj;
                  neigh.b = bj;
                  neigh_distances += deltaPrimeIJ(color, neigh, new_colors.data(), delta_prime_normalizer);
               }
            }

            bPrime.push_back(neigh_distances);
            // cerr << neigh_distances << endl;
            F[row][col] = neigh_distances;

            // cerr << neigh_distances << endl;
         }
      }
      cerr << "[DEBUG] b' vector calculated" << endl;

      //solve poisson
      double h1=2.0, h2=2.0, a1=1.0, a2=1.0;
      pde::types::boundary bdtype = pde::types::Neumann;
      double bdvalue=1.0;

      if(bdtype==pde::types::Neumann){
         bdvalue=pde::neumann_compat(F, a1, a2, h1, h2);
      }

      cerr << "[DEBUG] starting poisson solver" << endl;
      // pde::laplace(U,F,a1,a2,h1,h2,bdvalue,bdtype);
      double poi_trunc = pde::poisolve(U, F, a1, a2, h1, h2, bdvalue, bdtype, false);
      cerr << "[DEBUG] trunc=" << poi_trunc << endl;

      cerr << "[DEBUG] done poisson solver" << endl;

      // Normalize the result
      
      for (int i = 0; i < pSrc->GetHeight(); i++)
	   {
         for (int j = 0; j < pSrc->GetWidth(); j++)
         {
            if (U[i][j] < U_min)
               U_min = U[i][j];
            if (U[i][j] > U_max)
               U_max = U[i][j];
         }
      }

      // for (int i = 0; i < pSrc->GetHeight(); i++)
	   // {
      //    for (int j = 0; j < pSrc->GetWidth(); j++)
      //    {
      //       // U[i][j] = (U[i][j] - min) /(max - min);
      //    }
      // }
      
   }

   double origMeanL = getMeanL(pSourceData, pSrc->GetWidth(), pSrc->GetHeight());
   for (int i = 0; i < pSrc->GetHeight(); i++)
	{
		for (int j = 0; j < pSrc->GetWidth(); j++)
		{
         // cerr << U[i][j] << " ";
			// *pDestinationData++ = *pSourceData++;
			// *pDestinationData++ = *pSourceData++;
			// *pDestinationData++ = *pSourceData++;

         if(dJustRecolor){
            t_Lab new_color = new_colors[getClosestColorLab(*pSourceData++, *pSourceData++, *pSourceData++, new_colors.data(), new_colors.size())];
            // cerr << "[DEBUG] color=" << new_color.r << "(" << "-" << ") " << new_color.g << "(" << "-" << ") " << new_color.b << "(" << "-" << ")" << endl;
            *pDestinationData++ = new_color.L;
            *pDestinationData++ = new_color.a;
            *pDestinationData++ = new_color.b;
         }else{
            // *pDestinationData++ = U[i][j];
            // *pDestinationData++ = 0.0;
            // *pDestinationData++ = 0.0;

            double L;
            if(dNormalize == 1){
               // cerr << "[DEBUG] normalizing" << endl;
               L = (U[i][j] - U_min) / (U_max - U_min);
            }else if(dNormalize == 2){
               // cerr << "[DEBUG] offseting by mean L" << endl;
               L = U[i][j] + origMeanL;
            }else{
               // cerr << "[DEBUG] no normalization" << endl;
               L = U[i][j];
            }
            // L = U[i][j] / (U_max - U_min);
            setDataByCoords(i, j, pDestinationData, pSrc->GetWidth(), L, 0.0, 0.0);
         }
		}
      // cerr << endl;
	}

   cerr << "[DEBUG] Done" << endl;

	pDst->Convert(TMO_RGB);
	return 0;
}

// /* --------------------------------------------------------------------------- *
//  * This overloaded function is an implementation of your tone mapping operator *
//  * --------------------------------------------------------------------------- */
// int TMOxkocal00::Transform()
// {
//    cout << "[DEBUG] Starting" << endl;
// 	// Source image is stored in local parameter pSrc
// 	// Destination image is in pDst

// 	// Initialy images are in RGB format, but you can
// 	// convert it into other format
// 	pSrc->Convert(TMO_Yxy); // This is format of Y as luminance
// 	pDst->Convert(TMO_Yxy); // x, y as color information

// 	double *pSourceData = pSrc->GetData();		// You can work at low level data
// 	double *pDestinationData = pDst->GetData(); // Data are stored in form of array
// 												// of three doubles representing
// 												// three colour components
// 	double pY, px, py;

// 	int i = 0;
// 	for (i = 0; i < pSrc->GetHeight(); i++)
// 	{
// 		pSrc->ProgressBar(i, pSrc->GetHeight()); // You can provide progress bar
// 		for (int j = 0; j < pSrc->GetWidth(); j++)
// 		{
// 			// get data
//          getDataByCoords(i, j, pSourceData, pSrc->GetWidth(), &pY, &px, &py);

//          // average luminance in an 8 neighbourhood
//          int neighSide = (int) std::sqrt((float)dNeighbourhood);
//          double avgY = 0.0;
//          int missed_neighbours = 0;
//          for (int ii = -neighSide; ii <= neighSide; ii++){
//             for (int jj = -neighSide; jj <= neighSide; jj++){
//                // if (ii == 0 && jj == 0){
//                //    continue;
//                // }
//                if (i + ii < 0 || i + ii >= pSrc->GetHeight() || j+jj < 0 || j+jj >= pSrc->GetWidth()){
//                   missed_neighbours++;
//                   continue;
//                }

//                double tmpY, tmpx, tmpy;
//                getDataByCoords(i+ii, j+jj, pSourceData, pSrc->GetWidth(), &tmpY, &tmpx, &tmpy);

//                avgY += tmpY;
//                // px += tmpx;
//                // py += tmpy;
//             }
//          }

//          float neigh_count = (float) (neighSide - missed_neighbours);
//          pY = avgY / neigh_count;
//          // px /= neigh_count;
//          // py /= neigh_count;

// 			// and store results to the destination image
//          setDataByCoords(i, j, pDestinationData, pSrc->GetWidth(), pY, px, py);
// 		}
// 	}
// 	pSrc->ProgressBar(i, pSrc->GetHeight());
// 	pDst->Convert(TMO_RGB);
// 	return 0;
// }

void TMOxkocal00::getPixelByOffset(int y, int x, int dy, int dx, double *pixel){
   
}

