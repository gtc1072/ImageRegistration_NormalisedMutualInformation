#ifndef CLASS_NMI_REGISTRATION_BASEON_PARZON_WINDOW_HHHHHHH_____
#define CLASS_NMI_REGISTRATION_BASEON_PARZON_WINDOW_HHHHHHH_____

#include <vector>
#include <random>

#define R_PI 3.141592653589793238462643383279

enum TransformType
{
	translation = 0,
	rigid,
	affine,
	perspective,
	freeform
};

enum SamplerType
{
	random = 0,
	uniform
};

enum InterpolationType
{
	nn = 0,
	bilinear,
	cubic,
	bsplie
};

enum ParzonWindowKernelType
{
	gauss = 0,
	bspline,
	polynominal
};

struct ROIRegion
{
	int x, y, width, height;
};

class CnmiRegistration
{
public:
	CnmiRegistration();
	~CnmiRegistration();
	bool initRegistration(float* pFix, float* pMoved, int width, int height, SamplerType sampleType, TransformType transformType, InterpolationType interpolationType, ROIRegion roiFix, ROIRegion roiMove, int samplerCount,
		int iterCount, double iterEps, int bins);
	bool initRegistration(float* pFix, float* pMoved, int width, int height, SamplerType sampleType, TransformType transformType, InterpolationType interpolationType, unsigned char* pFixMask, unsigned char* pMoveMask, int samplerCount,
		int iterCount, double iterEps, int bins);
	std::vector<std::vector<float>> getTransformMatrix();
	bool setTransformMatrix(std::vector<std::vector<float>> matrix);
	void setLearnStepandMoment(double step, double moment);
	bool runRegistration();
private:
	void getMovedImageGradient();
	template<typename T>
	void freeImage(T *ptrImg);
	std::vector<std::pair<int, int>> RandomSampler(int count);
	std::vector<std::pair<int, int>> UniformSampler(int count);
	void initSampler();
	void initTransformMatirx();
	void initMask(ROIRegion roiFix, ROIRegion roiMove);
	void initIntensityRange();
	bool runRegistration_randomSampler();
	bool runRegistration_uniformSampler();
	bool updateSampler();
	std::pair<float, float> transformPoint(std::pair<int, int> in);
	std::vector<double> getFixImageProbabilityEstimation();
	std::vector<double> getMoveImageProbabilityEstimation();
	std::vector<std::vector<double>> getJointProbabilityEstimation(std::vector<double> &histFix, std::vector<double> &histMove);
	void getMoveImageminmax(int top, int bottom, int left, int right, float* pmin, float* pmax);
	float getNNInterpolationValue(float x, float y);
	float getBilinearInterpolationValue(float x, float y); 
	float getCubicInterpolationValue(float x, float y);
	float getBsplieInterpolationValue(float x, float y);
	std::pair<float, float> getNNInterpolationValueOfGradient(float x, float y);
	std::pair<float, float> getBilinearInterpolationValueOfGradient(float x, float y);
	std::pair<float, float> getCubicInterpolationValueOfGradient(float x, float y);
	std::pair<float, float> getBsplieInterpolationValueOfGradient(float x, float y);
	void getComputeResult(std::vector<double> &histFix, std::vector<double> &histMove, std::vector<std::vector<double>> &histJoint, double* B, double* A);
	void getComputeDerivate(std::vector<double> &histFix, std::vector<double> &histMove, std::vector<std::vector<double>> &histJoint, double B, double A);
	void updateTransform();
private:
	std::vector<std::vector<float>> m_transformMatrix;
	std::vector<std::vector<float>> m_transformMatrixDerivates;
	std::vector<std::vector<float>> m_preTransformMatrixDerivates;
	std::vector<std::vector<float>> m_rigidTransformMatrix;
	SamplerType m_sampleType;
	TransformType m_transformType; 
	InterpolationType m_interpolationType;
	int m_width, m_height, m_sampleCount;
	float* m_ptrFix;
	float* m_ptrMoved;
	float* m_ptrGradientX;
	float* m_ptrGradientY;
	ROIRegion m_roiFixRegion, m_roiMoveRegion;
	//ROIRegion m_roiRegion;
	unsigned char* m_ptrFixMask;
	unsigned char* m_ptrMoveMask;
	int m_fixBinSize;
	int m_moveBinSize;
	int m_fixMinBinIndex;
	int m_moveMinBinIndex;
	bool m_bFixUseMask;
	bool m_bMoveUseMask;
	bool m_bFixMaskSquare;
	bool m_bMoveMaskSquare;

	float m_fixMinValue;
	float m_fixMaxValue;
	float m_fixBinWidth;
	float m_fixRealMin;

	float m_moveMinValue;
	float m_moveMaxValue;
	float m_moveBinWidth;
	float m_moveRealMin;

	std::vector<std::pair<int, int>> m_samplePoint;
	std::vector<std::pair<int, int>> m_tempSamplePoint;
	std::vector<std::pair<float, float>> m_tempMovedSamplePoint;
	int m_tempSampleCount;
	int m_iterCount;
	double m_iterEps;
	int m_binsCount;
	double m_learnStep;
	double m_learnMoment;
	double* m_ptrSampletoHistValue;
};

#endif