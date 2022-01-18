#include "stdafx.h"
#include "CnmiRegistration.h"
#include <algorithm>
#include <numeric>
#include <functional>

#define LN2 0.693147181

CnmiRegistration::CnmiRegistration()
{
	m_ptrGradientX = m_ptrGradientY = m_ptrFix = m_ptrMoved = nullptr;
	m_ptrFixMask = m_ptrMoveMask = nullptr;
	m_fixBinSize = m_moveBinSize = m_fixMinBinIndex = m_moveMinBinIndex;
	m_bFixUseMask = m_bMoveUseMask = false;
	m_learnStep = 1.0e-2;
	m_learnMoment = 0.5;
	m_ptrSampletoHistValue = nullptr;
}

CnmiRegistration::~CnmiRegistration()
{
	freeImage(m_ptrGradientX);
	freeImage(m_ptrGradientY);
	freeImage(m_ptrFixMask);
	freeImage(m_ptrMoveMask);
	freeImage(m_ptrSampletoHistValue);
}

bool CnmiRegistration::initRegistration(float* pFix, float* pMoved, int width, int height, SamplerType sampleType, TransformType transformType, InterpolationType interpolationType, ROIRegion roiFix, ROIRegion roiMove, int samplerCount,
	int iterCount, double iterEps, int bins)
{
	if (pFix && pMoved && width > 100 && height > 100 && roiFix.x >= 0 && roiFix.y >= 0 && (roiFix.x + roiFix.width) <= width && (roiFix.y + roiFix.height) <= height && roiFix.width > 5 && roiFix.height > 5
		&& roiMove.x >= 0 && roiMove.y >= 0 && (roiMove.x + roiMove.width) <= width && (roiMove.y + roiMove.height) <= height && roiMove.width > 5 && roiMove.height > 5
		&& samplerCount > 100 && iterCount > 50 && iterEps > 0 && bins > 0)
	{
		freeImage(m_ptrGradientX);
		freeImage(m_ptrGradientY);
		freeImage(m_ptrFixMask);
		freeImage(m_ptrMoveMask);
		freeImage(m_ptrSampletoHistValue);
		m_ptrGradientX = new float[width * height];
		memset(m_ptrGradientX, 0x00, width * height * sizeof(float));
		m_ptrGradientY = new float[width * height];
		memset(m_ptrGradientY, 0x00, width * height * sizeof(float));
		m_ptrSampletoHistValue = new double[samplerCount * (bins * 2 + 2)];
		m_ptrFix = pFix;
		m_ptrMoved = pMoved;
		m_width = width;
		m_height = height;
		m_sampleType = sampleType;
		m_transformType = transformType;
		m_interpolationType = interpolationType;
		m_bFixUseMask = true;
		m_bFixMaskSquare = true;
		m_roiFixRegion = roiFix;
		if (roiFix.x == 0 && roiFix.y == 0 && roiFix.width == width && roiFix.height == height)
		{
			m_bFixUseMask = false;
		}
		m_bMoveUseMask = true;
		m_bMoveMaskSquare = true;
		m_roiMoveRegion = roiMove;
		if (roiMove.x == 0 && roiMove.y == 0 && roiMove.width == width && roiMove.height == height)
		{
			m_bMoveUseMask = false;
		}
		m_sampleCount = samplerCount;
		m_iterCount = iterCount;
		m_iterEps = iterEps;
		m_binsCount = bins;
		if (m_sampleCount > (m_width * m_height) / 10)
		{
			m_sampleCount = (m_width * m_height) / 10;
		}
		initMask(roiFix, roiMove);
		initIntensityRange();
		initTransformMatirx();
		getMovedImageGradient();
		return true;
	}
	return false;
}

bool CnmiRegistration::initRegistration(float* pFix, float* pMoved, int width, int height, SamplerType sampleType, TransformType transformType, InterpolationType interpolationType, unsigned char* pFixMask, unsigned char* pMoveMask, int samplerCount,
	int iterCount, double iterEps, int bins)
{
	if (pFix && pMoved && width > 100 && height > 100 && samplerCount > 100 && iterCount > 50 && iterEps > 0 && bins > 0)
	{
		freeImage(m_ptrGradientX);
		freeImage(m_ptrGradientY);
		freeImage(m_ptrFixMask);
		freeImage(m_ptrMoveMask);
		freeImage(m_ptrSampletoHistValue);
		m_ptrGradientX = new float[width * height];
		memset(m_ptrGradientX, 0x00, width * height * sizeof(float));
		m_ptrGradientY = new float[width * height];
		memset(m_ptrGradientY, 0x00, width * height * sizeof(float));
		m_ptrSampletoHistValue = new double[samplerCount * (bins * 2 + 2)];
		m_ptrFix = pFix;
		m_ptrMoved = pMoved;
		m_width = width;
		m_height = height;
		m_sampleType = sampleType;
		m_transformType = transformType;
		m_interpolationType = interpolationType;
		m_bFixUseMask = false;
		m_bFixMaskSquare = false;
		if (pFixMask)
		{
			m_bFixUseMask = true;
		}
		m_bMoveUseMask = false;
		m_bMoveMaskSquare = false;
		if (pMoveMask)
		{
			m_bMoveUseMask = true;
		}
		m_sampleCount = samplerCount;
		m_iterCount = iterCount;
		m_iterEps = iterEps;
		m_binsCount = bins;
		if (m_sampleCount > (m_width * m_height) / 10)
		{
			m_sampleCount = (m_width * m_height) / 10;
		}
		//initMask(roiFix, roiMove);
		m_ptrFixMask = pFixMask;
		m_ptrMoveMask = pMoveMask;
		initTransformMatirx();
		initIntensityRange();
		getMovedImageGradient();
		return true;
	}
	return false;
}

bool CnmiRegistration::runRegistration()
{
	if (m_sampleType == random)
		return runRegistration_randomSampler();
	return runRegistration_uniformSampler();
}

void CnmiRegistration::getComputeResult(std::vector<double> &histFix, std::vector<double> &histMove, std::vector<std::vector<double>> &histJoint, double* B, double* A)
{
	*B = *A = 0.0;
	for (int j = 0; j < m_binsCount; ++j)
	{
		*B += (-histFix[j] * std::log2(histFix[j]) - histMove[j] * std::log2(histMove[j]));
	}
	for (int j = 0; j < m_binsCount; ++j)
	{
		for (int i = 0; i < m_binsCount; ++i)
		{
			*A += (-histJoint[j][i] * std::log2(histJoint[j][i]));
		}
	}
}

void CnmiRegistration::getComputeDerivate(std::vector<double> &histFix, std::vector<double> &histMove, std::vector<std::vector<double>> &histJoint, double B, double A)
{
	if (m_transformType == translation)
	{
		m_transformMatrixDerivates[0][0] = m_transformMatrixDerivates[0][1] = 0.0;
		//m_preTransformMatrixDerivates[0][0] = m_preTransformMatrixDerivates[0][1] = 0.0;
		const double p = sqrt(2.0 * R_PI);
		double h = 0.4;
		double h2 = h*h;
		double h3 = h2*h;
		std::function<float(float, float)> fun;
		std::function<std::pair<float, float>(float, float)> fun_g;
		if (m_interpolationType == nn)
		{
			fun = std::bind(&CnmiRegistration::getNNInterpolationValue, this, std::placeholders::_1, std::placeholders::_2);
			fun_g = std::bind(&CnmiRegistration::getNNInterpolationValueOfGradient, this, std::placeholders::_1, std::placeholders::_2);
		}
		else if (m_interpolationType == bilinear)
		{
			fun = std::bind(&CnmiRegistration::getBilinearInterpolationValue, this, std::placeholders::_1, std::placeholders::_2);
			fun_g = std::bind(&CnmiRegistration::getBilinearInterpolationValueOfGradient, this, std::placeholders::_1, std::placeholders::_2);
		}
		else if (m_interpolationType == cubic)
		{
			fun = std::bind(&CnmiRegistration::getCubicInterpolationValue, this, std::placeholders::_1, std::placeholders::_2);
			fun_g = std::bind(&CnmiRegistration::getCubicInterpolationValueOfGradient, this, std::placeholders::_1, std::placeholders::_2);
		}
		else if (m_interpolationType == bsplie)
		{
			fun = std::bind(&CnmiRegistration::getBsplieInterpolationValue, this, std::placeholders::_1, std::placeholders::_2);
			fun_g = std::bind(&CnmiRegistration::getBsplieInterpolationValueOfGradient, this, std::placeholders::_1, std::placeholders::_2);
		}
		int sampletoHistLength = m_binsCount * 2 + 2;
		for (int j = 0; j < m_binsCount; ++j)
		{
			for (int i = 0; i < m_binsCount; ++i)
			{
				double partion = A * std::log2(histMove[i]) - B * std::log2(histJoint[j][i]);
				partion *= (1.0 / m_tempSampleCount) /*/ m_moveBinWidth*/;
				//m_transformMatrixDerivates[0][0] += 
				double sum_x = 0.0;
				double sum_y = 0.0;
				for (int k = 0; k < m_tempSampleCount; ++k)
				{
					double *pSampleToHist = m_ptrSampletoHistValue + k * sampletoHistLength;
					float value_f = pSampleToHist[0];
					float value_m = pSampleToHist[1];
					//int x_f = m_tempSamplePoint[k].first;
					//int y_f = m_tempSamplePoint[k].second;
					//float value_f = m_ptrFix[y_f*m_width + x_f];
					//value_f = value_f / m_fixBinWidth - m_fixRealMin;
					//int val_f = (int)value_f;
					//if (val_f < 2.0) val_f = 2.0;
					//if (val_f > m_binsCount - 3) val_f = m_binsCount - 3;
					if (abs(value_f - j) > 2) continue;
					//value_f = val_f;
					float x_m = m_tempMovedSamplePoint[k].first;
					float y_m = m_tempMovedSamplePoint[k].second;
					//float value_m = fun(x_m, y_m);
					//value_m = value_m / m_moveBinWidth - m_moveRealMin;
					//int val_m = (int)value_m;
					//if (val_m < 2.0) val_m = 2.0;
					//if (val_m > m_binsCount - 3) val_m = m_binsCount - 3;
					if (abs(value_m - i) > 2) continue;
					//value_m = val_m;
					double temp_value_x = 0.0;
					double temp_value_y = 0.0;
					//temp_value_x = 1.0 / (p * h * 1.084997) * exp(-0.5 * (j - value_f)*(j - value_f) / h2);
					//temp_value_x *= -1.0 * (i - value_m) / (p * h3 * 1.084997) * exp(-0.5 * (i - value_m) * (i - value_m) / h2);
					temp_value_y = temp_value_x = pSampleToHist[j + 2] * pSampleToHist[i + m_binsCount + 2];
					std::pair<float, float> gxy = fun_g(x_m, y_m);
					temp_value_x *= (-1.0 * gxy.first);
					temp_value_y *= (-1.0 * gxy.second);
					sum_x += temp_value_x;
					sum_y += temp_value_y;
				}
				m_transformMatrixDerivates[0][0] += (sum_x * partion);
				m_transformMatrixDerivates[0][1] += (sum_y * partion);
			}
		}
		m_transformMatrixDerivates[0][0] /= (-1.0 * A*A);
		m_transformMatrixDerivates[0][1] /= (-1.0 * A*A);
	}
	else if (m_transformType == rigid)
	{
		m_transformMatrixDerivates[0][0] = m_transformMatrixDerivates[0][1] = m_transformMatrixDerivates[0][2] = m_transformMatrixDerivates[0][3] = m_transformMatrixDerivates[0][4] = 0.0;
		const double p = sqrt(2.0 * R_PI);
		double h = 0.4;
		double h2 = h*h;
		double h3 = h2*h;
		std::function<float(float, float)> fun;
		std::function<std::pair<float, float>(float, float)> fun_g;
		if (m_interpolationType == nn)
		{
			fun = std::bind(&CnmiRegistration::getNNInterpolationValue, this, std::placeholders::_1, std::placeholders::_2);
			fun_g = std::bind(&CnmiRegistration::getNNInterpolationValueOfGradient, this, std::placeholders::_1, std::placeholders::_2);
		}
		else if (m_interpolationType == bilinear)
		{
			fun = std::bind(&CnmiRegistration::getBilinearInterpolationValue, this, std::placeholders::_1, std::placeholders::_2);
			fun_g = std::bind(&CnmiRegistration::getBilinearInterpolationValueOfGradient, this, std::placeholders::_1, std::placeholders::_2);
		}
		else if (m_interpolationType == cubic)
		{
			fun = std::bind(&CnmiRegistration::getCubicInterpolationValue, this, std::placeholders::_1, std::placeholders::_2);
			fun_g = std::bind(&CnmiRegistration::getCubicInterpolationValueOfGradient, this, std::placeholders::_1, std::placeholders::_2);
		}
		else if (m_interpolationType == bsplie)
		{
			fun = std::bind(&CnmiRegistration::getBsplieInterpolationValue, this, std::placeholders::_1, std::placeholders::_2);
			fun_g = std::bind(&CnmiRegistration::getBsplieInterpolationValueOfGradient, this, std::placeholders::_1, std::placeholders::_2);
		}
		int sampletoHistLength = m_binsCount * 2 + 2;

		double cur_cx = m_transformMatrix[0][1];
		double cur_cy = m_transformMatrix[0][2];
		double cur_tx = m_transformMatrix[0][3];
		double cur_ty = m_transformMatrix[0][4];

		double cur_sin_theta = std::sin(m_transformMatrix[0][0]);
		double cur_cos_theta = std::cos(m_transformMatrix[0][0]);

		for (int j = 0; j < m_binsCount; ++j)
		{
			for (int i = 0; i < m_binsCount; ++i)
			{
				double partion = A * std::log2(histMove[i]) - B * std::log2(histJoint[j][i]);
				partion *= (1.0 / m_tempSampleCount);
				double sum_theta = 0.0;
				double sum_cx = 0.0;
				double sum_cy = 0.0;
				double sum_tx = 0.0;
				double sum_ty = 0.0;
				for (int k = 0; k < m_tempSampleCount; ++k)
				{
					double *pSampleToHist = m_ptrSampletoHistValue + k * sampletoHistLength;
					float value_f = pSampleToHist[0];
					float value_m = pSampleToHist[1];
					if (abs(value_f - j) > 2) continue;
					float x_m = m_tempMovedSamplePoint[k].first;
					float y_m = m_tempMovedSamplePoint[k].second;
					int x_f = m_tempSamplePoint[k].first;
					int y_f = m_tempSamplePoint[k].second;
					if (abs(value_m - i) > 2) continue;
					double temp_value_theta = 0.0;
					double temp_value_cx = 0.0;
					double temp_value_cy = 0.0;
					double temp_value_tx = 0.0;
					double temp_value_ty = 0.0;

					temp_value_theta = temp_value_cx = temp_value_cy = temp_value_tx = temp_value_ty = pSampleToHist[j + 2] * pSampleToHist[i + m_binsCount + 2];
					std::pair<float, float> gxy = fun_g(x_m, y_m);


					temp_value_theta *= gxy.first * ((x_f - cur_cx) * cur_sin_theta + (y_f - cur_cy) * cur_cos_theta) + gxy.second * ((y_f - cur_cy) * cur_sin_theta - (x_f - cur_cx) * cur_cos_theta);
					temp_value_cx *= gxy.first * cur_cos_theta + gxy.second * cur_sin_theta;
					temp_value_cy *= (-1.0 * gxy.first) * cur_sin_theta + gxy.second * cur_cos_theta;
					temp_value_tx *= (-1.0 * gxy.first);
					temp_value_ty *= (-1.0 * gxy.second);
					sum_theta += temp_value_theta;
					sum_cx += temp_value_cx;
					sum_cy += temp_value_cy;
					sum_tx += temp_value_tx;
					sum_ty += temp_value_ty;
				}
				m_transformMatrixDerivates[0][0] += (sum_theta * partion);
				m_transformMatrixDerivates[0][1] += (sum_cx * partion);
				m_transformMatrixDerivates[0][2] += (sum_cy * partion);
				m_transformMatrixDerivates[0][3] += (sum_tx * partion);
				m_transformMatrixDerivates[0][4] += (sum_ty * partion);
			}
		}
		m_transformMatrixDerivates[0][0] /= (-1.0 * A*A);
		m_transformMatrixDerivates[0][1] /= (-1.0 * A*A);
		m_transformMatrixDerivates[0][2] /= (-1.0 * A*A);
		m_transformMatrixDerivates[0][3] /= (-1.0 * A*A);
		m_transformMatrixDerivates[0][4] /= (-1.0 * A*A);
	}
	else if (m_transformType == affine)
	{

	}
	else if (m_transformType == perspective)
	{

	}
	else if (m_transformType == freeform)
	{

	}
}

void CnmiRegistration::setLearnStepandMoment(double step, double moment)
{
	if (step > 0)
		m_learnStep = step;
	if (moment > 0 && moment < 1.0)
		m_learnMoment = moment;
}

void CnmiRegistration::updateTransform()
{
	if (m_transformType == translation)
	{
		m_preTransformMatrixDerivates[0][0] = m_preTransformMatrixDerivates[0][0] * m_learnMoment + m_transformMatrixDerivates[0][0];
		m_preTransformMatrixDerivates[0][1] = m_preTransformMatrixDerivates[0][1] * m_learnMoment + m_transformMatrixDerivates[0][1];
		m_transformMatrix[0][0] += m_preTransformMatrixDerivates[0][0] * m_learnStep;
		m_transformMatrix[0][1] += m_preTransformMatrixDerivates[0][1] * m_learnStep;
		//printf("%.2f, %.2f\n", m_transformMatrix[0][0], m_transformMatrix[0][1]);
	}
	else if (m_transformType == rigid)
	{
		m_preTransformMatrixDerivates[0][0] = m_preTransformMatrixDerivates[0][0] * m_learnMoment + m_transformMatrixDerivates[0][0];
		m_preTransformMatrixDerivates[0][1] = m_preTransformMatrixDerivates[0][1] * m_learnMoment + m_transformMatrixDerivates[0][1];
		m_preTransformMatrixDerivates[0][2] = m_preTransformMatrixDerivates[0][2] * m_learnMoment + m_transformMatrixDerivates[0][2];
		m_preTransformMatrixDerivates[0][3] = m_preTransformMatrixDerivates[0][3] * m_learnMoment + m_transformMatrixDerivates[0][3];
		m_preTransformMatrixDerivates[0][4] = m_preTransformMatrixDerivates[0][4] * m_learnMoment + m_transformMatrixDerivates[0][4];
		m_transformMatrix[0][0] += m_preTransformMatrixDerivates[0][0] * m_learnStep * 0.001;
		m_transformMatrix[0][1] += m_preTransformMatrixDerivates[0][1] * m_learnStep;
		m_transformMatrix[0][2] += m_preTransformMatrixDerivates[0][2] * m_learnStep;
		m_transformMatrix[0][3] += m_preTransformMatrixDerivates[0][3] * m_learnStep;
		m_transformMatrix[0][4] += m_preTransformMatrixDerivates[0][4] * m_learnStep;
		printf("angle: %.6f, c(%.2f, %.2f) t(%.2f, %.2f)\n", m_transformMatrix[0][0], m_transformMatrix[0][1], m_transformMatrix[0][2], m_transformMatrix[0][3], m_transformMatrix[0][4]);
		double radian = m_transformMatrix[0][0];
		double cx = m_transformMatrix[0][1];
		double cy = m_transformMatrix[0][2];
		double tx = m_transformMatrix[0][3];
		double ty = m_transformMatrix[0][4];
		m_rigidTransformMatrix[0][0] = m_rigidTransformMatrix[1][1] = std::cos(radian);
		m_rigidTransformMatrix[0][1] = -std::sin(radian);
		m_rigidTransformMatrix[1][0] = std::sin(radian);
		m_rigidTransformMatrix[2][2] = 1.0;
		m_rigidTransformMatrix[0][2] = -cx * std::cos(radian) + cy * std::sin(radian) + tx;
		m_rigidTransformMatrix[1][2] = -cx * std::sin(radian) - cy * std::cos(radian) + ty;

	}
	else if (m_transformType == affine)
	{

	}
	else if (m_transformType == perspective)
	{

	}
	else if (m_transformType == freeform)
	{

	}
}

bool CnmiRegistration::runRegistration_randomSampler()
{
	int loop = 0; 
	double err = m_iterEps + 1.0;
	double B, A;
	while (loop < m_iterCount && err > m_iterEps)
	{
		initSampler();
		updateSampler();
		m_tempSampleCount = m_tempSamplePoint.size();
		if (m_tempSampleCount < (int)(m_sampleCount * 0.4))
		{
			printf("overlap region smaller than expected.\n");
			return false;
		}
		//std::vector<double> fixedProbability = getFixImageProbabilityEstimation();
		//std::vector<double> movedProbability = getMoveImageProbabilityEstimation();
		memset(m_ptrSampletoHistValue, 0x0, m_sampleCount * (m_binsCount * 2 + 2) * sizeof(double));
		std::vector<double> histFix, histMove;
		std::vector<std::vector<double>> jointProbability = getJointProbabilityEstimation(histFix, histMove);
		B = A = 0.0;
		getComputeResult(histFix, histMove, jointProbability, &B, &A);
		getComputeDerivate(histFix, histMove, jointProbability, B, A);
		std::vector<std::vector<float>> preTransform = m_transformMatrix;
		//printf("NMI=%.3f\n", B / A);
		updateTransform();
		err = 0.0;
		for (int j = 0; j < preTransform.size(); ++j)
		{
			for (int i = 0; i < preTransform[j].size(); ++i)
			{
				err += abs(preTransform[j][i] - m_transformMatrix[j][i]);
			}
		}
		loop++;
	}
	return true;
}

bool CnmiRegistration::runRegistration_uniformSampler()
{
	return false;
}

std::vector<std::vector<float>> CnmiRegistration::getTransformMatrix()
{
	switch (m_transformType)
	{
	case translation:
		return m_transformMatrix;
	case rigid:
		return m_rigidTransformMatrix;
	case affine:
		break;
	case perspective:
		break;
	case freeform:
		break;
	default:
		break;
	}
	return m_transformMatrix;
}

void CnmiRegistration::initIntensityRange()
{
	m_fixMinValue = 1.0e+20;
	m_fixMaxValue = -1.0e+20;
	m_moveMinValue = 1.0e+20;
	m_moveMaxValue = -1.0e+20;
	const int padding = 2;
	if (m_bFixUseMask)
	{
		if (m_bFixMaskSquare)
		{
			int top = m_roiFixRegion.y;
			int bottom = m_roiFixRegion.y + m_roiFixRegion.height;
			int left = m_roiFixRegion.x;
			int right = m_roiFixRegion.x + m_roiFixRegion.width;
			for (int j = top; j < bottom; ++j)
			{
				for (int i = left; i < right; ++i)
				{
					float v = m_ptrFix[j*m_width + i];
					if (v < m_fixMinValue) m_fixMinValue = v;
					if (v > m_fixMaxValue) m_fixMaxValue = v;
				}
			}
		}
		else
		{
			int step = (m_width + 3) / 4 * 4;
			for (int j = 0; j < m_height; ++j)
			{
				for (int i = 0; i < m_width; ++i)
				{
					if (m_ptrFixMask[j * step + i] == 255)
					{
						float v = m_ptrFix[j*m_width + i];
						if (v < m_fixMinValue) m_fixMinValue = v;
						if (v > m_fixMaxValue) m_fixMaxValue = v;
					}
				}
			}
		}
	}
	else
	{
		for (int j = 0; j < m_height; ++j)
		{
			for (int i = 0; i < m_width; ++i)
			{
				float v = m_ptrFix[j*m_width + i];
				if (v < m_fixMinValue) m_fixMinValue = v;
				if (v > m_fixMaxValue) m_fixMaxValue = v;
			}
		}
	}
	
	m_fixBinWidth = (m_fixMaxValue - m_fixMinValue) / ((float)(m_binsCount - 2 * padding));
	m_fixRealMin = m_fixMinValue / m_fixBinWidth - padding;

	if (m_bMoveUseMask)
	{
		if (m_bMoveMaskSquare)
		{
			int top = m_roiMoveRegion.y;
			int bottom = m_roiMoveRegion.y + m_roiMoveRegion.height;
			int left = m_roiMoveRegion.x;
			int right = m_roiMoveRegion.x + m_roiMoveRegion.width;
			for (int j = top; j < bottom; ++j)
			{
				for (int i = left; i < right; ++i)
				{
					float v = m_ptrFix[j*m_width + i];
					if (v < m_moveMinValue) m_moveMinValue = v;
					if (v > m_moveMaxValue) m_moveMaxValue = v;
				}
			}
		}
		else
		{
			int step = (m_width + 3) / 4 * 4;
			for (int j = 0; j < m_height; ++j)
			{
				for (int i = 0; i < m_width; ++i)
				{
					if (m_ptrMoveMask[j * step + i] == 255)
					{
						float v = m_ptrMoved[j*m_width + i];
						if (v < m_moveMinValue) m_moveMinValue = v;
						if (v > m_moveMaxValue) m_moveMaxValue = v;
					}
				}
			}
		}
	}
	else
	{
		for (int j = 0; j < m_height; ++j)
		{
			for (int i = 0; i < m_width; ++i)
			{
				float v = m_ptrMoved[j*m_width + i];
				if (v < m_moveMinValue) m_moveMinValue = v;
				if (v > m_moveMaxValue) m_moveMaxValue = v;
			}
		}
	}

	m_moveBinWidth = (m_moveMaxValue - m_moveMinValue) / ((float)(m_binsCount - 2 * padding));
	m_moveRealMin = m_moveMinValue / m_moveBinWidth - padding;
}

void CnmiRegistration::initMask(ROIRegion roiFix, ROIRegion roiMove)
{
	freeImage(m_ptrFixMask);
	freeImage(m_ptrMoveMask);
	int step = (m_width + 3) / 4 * 4;
	if (m_bFixUseMask)
	{
		m_ptrFixMask = new unsigned char[step * m_height];
		memset(m_ptrFixMask, 0x00, step * m_height);
		for (int j = roiFix.y; j < roiFix.y + roiFix.height; ++j)
		{
			for (int i = roiFix.x; i < roiFix.x + roiFix.width; ++i)
			{
				m_ptrFixMask[j*step + i] = 255;
			}
		}
	}
	if (m_bMoveUseMask)
	{
		m_ptrMoveMask = new unsigned char[step * m_height];
		memset(m_ptrMoveMask, 0x00, step * m_height);
		for (int j = roiMove.y; j < roiMove.y + roiMove.height; ++j)
		{
			for (int i = roiMove.x; i < roiMove.x + roiMove.width; ++i)
			{
				m_ptrMoveMask[j*step + i] = 255;
			}
		}
	}
}

void CnmiRegistration::initSampler()
{
	if (m_sampleType == random)
	{
		m_samplePoint = RandomSampler(m_sampleCount);
	}
	else
	{
		m_samplePoint = UniformSampler(m_sampleCount);
		m_sampleCount = m_samplePoint.size();
	}
}

bool CnmiRegistration::updateSampler()
{
	bool isAllinside = true;
	m_tempSamplePoint.clear();
	m_tempMovedSamplePoint.clear();
	if (m_transformType == translation)
	{
		std::pair<float, float> out;
		if (m_bMoveUseMask)
		{
			if (m_bMoveMaskSquare)
			{
				for (std::vector<std::pair<int, int>>::iterator it = m_samplePoint.begin(); it != m_samplePoint.end(); ++it)
				{
					out.first = it->first + m_transformMatrix[0][0];
					out.second = it->second + m_transformMatrix[0][1];
					if (out.first < m_roiMoveRegion.x + 2 || out.first >= (m_roiMoveRegion.x + m_roiMoveRegion.width) - 3 || out.second < m_roiMoveRegion.y + 2 || out.second >= (m_roiMoveRegion.y + m_roiMoveRegion.height) - 3)
					{
						isAllinside = false;
						continue;
					}
					m_tempSamplePoint.push_back(*it);
					m_tempMovedSamplePoint.push_back(out);
				}
			}
			else
			{
				int step = (m_width + 3) / 4 * 4;
				for (std::vector<std::pair<int, int>>::iterator it = m_samplePoint.begin(); it != m_samplePoint.end(); ++it)
				{
					out.first = it->first + m_transformMatrix[0][0];
					out.second = it->second + m_transformMatrix[0][1];
					if (m_ptrMoveMask[std::lround(out.second) * step + std::lround(out.first)] == 0 || out.first < 2 || out.first >= m_width - 3 || out.second < 2 || out.second >= m_height - 3)
					{
						isAllinside = false;
						continue;
					}
					m_tempSamplePoint.push_back(*it);
					m_tempMovedSamplePoint.push_back(out);
				}
			}
		}
		else
		{
			for (std::vector<std::pair<int, int>>::iterator it = m_samplePoint.begin(); it != m_samplePoint.end(); ++it)
			{
				out.first = it->first + m_transformMatrix[0][0];
				out.second = it->second + m_transformMatrix[0][1];
				if (out.first < 2 || out.first >= m_width - 3 || out.second < 2 || out.second >= m_height - 3)
				{
					isAllinside = false;
					continue;
				}
				m_tempSamplePoint.push_back(*it);
				m_tempMovedSamplePoint.push_back(out);
			}
		}
	}
	else if (m_transformType == rigid)
	{
		std::pair<float, float> out;
		if (m_bMoveUseMask)
		{
			if (m_bMoveMaskSquare)
			{
				for (std::vector<std::pair<int, int>>::iterator it = m_samplePoint.begin(); it != m_samplePoint.end(); ++it)
				{
					out.first = it->first * m_rigidTransformMatrix[0][0] + it->second * m_rigidTransformMatrix[0][1] + m_rigidTransformMatrix[0][2];
					out.second = it->first * m_rigidTransformMatrix[1][0] + it->second * m_rigidTransformMatrix[1][1] + m_rigidTransformMatrix[1][2];
					if (out.first < m_roiMoveRegion.x + 2 || out.first >= (m_roiMoveRegion.x + m_roiMoveRegion.width) - 3 || out.second < m_roiMoveRegion.y + 2 || out.second >= (m_roiMoveRegion.y + m_roiMoveRegion.height) - 3)
					{
						isAllinside = false;
						continue;
					}
					m_tempSamplePoint.push_back(*it);
					m_tempMovedSamplePoint.push_back(out);
				}
			}
			else
			{
				int step = (m_width + 3) / 4 * 4;
				for (std::vector<std::pair<int, int>>::iterator it = m_samplePoint.begin(); it != m_samplePoint.end(); ++it)
				{
					out.first = it->first * m_rigidTransformMatrix[0][0] + it->second * m_rigidTransformMatrix[0][1] + m_rigidTransformMatrix[0][2];
					out.second = it->first * m_rigidTransformMatrix[1][0] + it->second * m_rigidTransformMatrix[1][1] + m_rigidTransformMatrix[1][2];
					if (m_ptrMoveMask[std::lround(out.second) * step + std::lround(out.first)] == 0 || out.first < 2 || out.first >= m_width - 3 || out.second < 2 || out.second >= m_height - 3)
					{
						isAllinside = false;
						continue;
					}
					m_tempSamplePoint.push_back(*it);
					m_tempMovedSamplePoint.push_back(out);
				}
			}
		}
		else
		{
			for (std::vector<std::pair<int, int>>::iterator it = m_samplePoint.begin(); it != m_samplePoint.end(); ++it)
			{
				out.first = it->first * m_rigidTransformMatrix[0][0] + it->second * m_rigidTransformMatrix[0][1] + m_rigidTransformMatrix[0][2];
				out.second = it->first * m_rigidTransformMatrix[1][0] + it->second * m_rigidTransformMatrix[1][1] + m_rigidTransformMatrix[1][2];
				if (out.first < 2 || out.first >= m_width - 3 || out.second < 2 || out.second >= m_height - 3)
				{
					isAllinside = false;
					continue;
				}
				m_tempSamplePoint.push_back(*it);
				m_tempMovedSamplePoint.push_back(out);
			}
		}
	}
	else if (m_transformType == affine)
	{

	}
	else if (m_transformType == perspective)
	{

	}
	else if (m_transformType == freeform)
	{

	}
	return isAllinside;
}

std::pair<float, float> CnmiRegistration::transformPoint(std::pair<int, int> in)
{
	std::pair<float, float> out;
	out.first = in.first + m_transformMatrix[0][0];
	out.second = in.second + m_transformMatrix[0][1];
	return out;
}

void CnmiRegistration::initTransformMatirx()
{
	m_transformMatrix.clear();
	m_transformMatrixDerivates.clear();
	m_preTransformMatrixDerivates.clear();
	std::vector<float> temp;
	switch (m_transformType)
	{
	case translation:
		temp.push_back(0.0);
		temp.push_back(0.0);
		m_transformMatrix.push_back(temp);
		m_transformMatrixDerivates.push_back(temp);
		m_preTransformMatrixDerivates.push_back(temp);
		break;
	case rigid:
		temp.push_back(0.0);//theta
		temp.push_back(m_width / 2);//cx
		temp.push_back(m_height / 2);//cy
		temp.push_back(m_width / 2);//tx
		temp.push_back(m_height / 2);//ty
		m_transformMatrix.push_back(temp);
		m_transformMatrixDerivates = std::vector<std::vector<float>>(1, std::vector<float>(5, 0));
		m_preTransformMatrixDerivates = std::vector<std::vector<float>>(1, std::vector<float>(5, 0));
		m_rigidTransformMatrix = std::vector<std::vector<float>>(3, std::vector<float>(3, 0));
		m_rigidTransformMatrix[0][0] = m_rigidTransformMatrix[1][1] = m_rigidTransformMatrix[2][2] = 1.0;
		//m_rigidTransformMatrix[0][2] = -m_width / 2;
		//m_rigidTransformMatrix[1][2] = -m_height / 2;
		break;
	case affine:
		break;
	case perspective:
		break;
	case freeform:
		break;
	default:
		break;
	}
}

template<typename T>
void CnmiRegistration::freeImage(T *ptrImg)
{
	if (ptrImg)
	{
		delete ptrImg;
		ptrImg = nullptr;
	}
}

void CnmiRegistration::getMovedImageGradient()
{
	double top = 0.0;
	double bottom = 0.0;
	double middle = 0.0;
	for (int j = 1; j < m_height - 1; ++j)
	{
		float* pTop = m_ptrMoved + (j - 1) * m_width;
		float* pMid = m_ptrMoved + j * m_width;
		float* pBottom = m_ptrMoved + (j + 1) * m_width;
		float* pX = m_ptrGradientX + j * m_width;
		float* pY = m_ptrGradientY + j * m_width;
		for (int i = 1; i < m_width - 1; ++i)
		{
			top = pTop[i + 1] - pTop[i - 1];
			middle = pMid[i + 1] - pMid[i - 1];
			bottom = pBottom[i + 1] - pBottom[i - 1];
			pX[i] = /*0.25 **/ (top + middle * 2.0 + bottom);
			top = pBottom[i - 1] - pTop[i - 1];
			middle = pBottom[i] - pTop[i];
			bottom = pBottom[i + 1] - pTop[i + 1];
			pY[i] = /*0.25 **/ (top + middle * 2.0 + bottom);
		}
	}
}

std::vector<std::pair<int, int>> CnmiRegistration::RandomSampler(int count)
{
	std::vector<std::pair<int, int>> samples;
	std::random_device rd;
	std::default_random_engine generator(rd());
	std::uniform_int_distribution<int> distribution_x(0, m_width - 1);
	std::uniform_int_distribution<int> distribution_y(0, m_height - 1);
	if (m_bFixUseMask)
	{
		if (m_bFixMaskSquare)
		{
			distribution_x = std::uniform_int_distribution<int>(m_roiFixRegion.x, m_roiFixRegion.x + m_roiFixRegion.width - 1);
			distribution_y = std::uniform_int_distribution<int>(m_roiFixRegion.y, m_roiFixRegion.y + m_roiFixRegion.height - 1);
			for (int i = 0; i < count; ++i)
			{
				samples.push_back(std::make_pair(distribution_x(generator), distribution_y(generator)));
			}
		}
		else
		{
			int step = (m_width + 3) / 4 * 4;
			for (int i = 0; i < count; ++i)
			{
				int x = distribution_x(generator);
				int y = distribution_y(generator);
				if (m_ptrFixMask[y * step + x] == 255)
					samples.push_back(std::make_pair(x, y));
			}
		}
	}
	else
	{
		for (int i = 0; i < count; ++i)
		{
			samples.push_back(std::make_pair(distribution_x(generator), distribution_y(generator)));
		}
	}
	return samples;
}

std::vector<std::pair<int, int>> CnmiRegistration::UniformSampler(int count)
{
	std::vector<std::pair<int, int>> samples;
	if (m_bFixUseMask)
	{
		if (m_bFixMaskSquare)
		{
			int deltaW = (int)std::floor(sqrt((double)(m_roiFixRegion.width * m_roiFixRegion.height) / (double)(count)));
			int wid_count = (int)std::floor((double)(m_roiFixRegion.width) / (double)(deltaW));
			int wid_count2 = (int)std::floor((double)(m_roiFixRegion.height) / (double)(deltaW));
			wid_count = std::min(wid_count, wid_count2);
			int offx = m_roiFixRegion.x + deltaW / 2;
			int offy = m_roiFixRegion.y + deltaW / 2;
			for (int i = 0; i < count; ++i)
			{
				int p_j = i / wid_count;
				int p_i = i % wid_count;
				p_j = p_j * deltaW + offy;
				p_i = p_i * deltaW + offx;
				if (p_j < m_roiFixRegion.y + m_roiFixRegion.height)
					samples.push_back(std::make_pair(p_i, p_j));
			}
		}
		else
		{
			//complex
		}
	}
	else
	{
		int deltaW = (int)std::floor(sqrt((double)(m_width*m_height) / (double)(count)));
		int wid_count = (int)std::floor((double)(m_width) / (double)(deltaW));
		int wid_count2 = (int)std::floor((double)(m_height) / (double)(deltaW));
		wid_count = std::min(wid_count, wid_count2);
		int offx = deltaW / 2;
		int offy = deltaW / 2;
		for (int i = 0; i < count; ++i)
		{
			int p_j = i / wid_count;
			int p_i = i % wid_count;
			p_j = p_j * deltaW + offy;
			p_i = p_i * deltaW + offx;
			if (p_j < m_height)
				samples.push_back(std::make_pair(p_i, p_j));
		}
	}
	return samples;
}

bool CnmiRegistration::setTransformMatrix(std::vector<std::vector<float>> matrix)
{
	if (m_transformType == translation)
	{
		if (matrix.empty())
			return false;
		if (matrix[0].size() < 2)
			return false;
		m_transformMatrix = matrix;
	}
	else if (m_transformType == rigid)
	{
		if (matrix.empty())
			return false;
		if (matrix[0].size() < 5)
			return false;
		m_rigidTransformMatrix = std::vector<std::vector<float>>(3, std::vector<float>(3, 0));
		double radian = matrix[0][0] / 180.0 * R_PI;
		double cx = matrix[0][1];
		double cy = matrix[0][2];
		double tx = matrix[0][3];
		double ty = matrix[0][4];
		m_rigidTransformMatrix[0][0] = m_rigidTransformMatrix[1][1] = std::cos(radian);
		m_rigidTransformMatrix[0][1] = -std::sin(radian);
		m_rigidTransformMatrix[1][0] = std::sin(radian);
		m_rigidTransformMatrix[2][2] = 1.0;
		m_rigidTransformMatrix[0][2] = -cx * std::cos(radian) + cy * std::sin(radian) + tx;
		m_rigidTransformMatrix[1][2] = -cx * std::sin(radian) - cy * std::cos(radian) + ty;
		m_transformMatrix = matrix;
		m_transformMatrix[0][0] = radian;
	}
	else if (m_transformType == affine)
	{

	}
	else if (m_transformType == perspective)
	{

	}
	else if (m_transformType == freeform)
	{

	}
	return true;
}

std::vector<double> CnmiRegistration::getFixImageProbabilityEstimation()
{
	std::vector<double> prob(m_binsCount, 0.0);
	const double p = sqrt(2.0 * R_PI);
	//double h = 0.73;
	double h = 0.4;
	////
	for (int i = 0; i < m_tempSampleCount; ++i)
	{
		int x = m_tempSamplePoint[i].first;
		int y = m_tempSamplePoint[i].second;
		float value = m_ptrFix[y*m_width + x];
		value = value / m_fixBinWidth - m_fixRealMin;
		int val = (int)value;
		if (val < 2.0) val = 2.0;
		if (val > m_binsCount - 3) val = m_binsCount - 3;
		value = val;
		for (int k = 0; k < m_binsCount; ++k)
		{
			//0.73
			//prob[k] += 1.0 / (m_tempSampleCount * p * h) * exp(-0.5 * (k - value)*(k - value) / (h * h));
			//0.4
			prob[k] += 1.0 / (m_tempSampleCount * p * h * 1.084997) * exp(-0.5 * (k - value)*(k - value) / (h * h));
		}
	}

	for (int k = 0; k < m_binsCount; ++k)
	{
		if (prob[k] < 1.0e-50)
		{
			prob[k] = 1.0e-50;
		}
	}

	return prob;
}

float CnmiRegistration::getNNInterpolationValue(float x, float y)
{
	float ret = m_ptrMoved[std::lroundf(y) * m_width + std::lroundf(x)];
	return ret;
}

float CnmiRegistration::getBilinearInterpolationValue(float x, float y)
{
	float ret = 0.0;
	int ix = lround(x);
	int iy = lround(y);
	ret = ((1.0 + ix - x) * m_ptrMoved[iy * m_width + ix] + (x - ix) * m_ptrMoved[iy * m_width + ix + 1]) * (1.0 + iy - y)
		+ ((1.0 + ix - x) * m_ptrMoved[(iy + 1) * m_width + ix] + (x - ix) * m_ptrMoved[(iy + 1) * m_width + ix + 1]) * (y - iy);
	return ret;
}

inline float cubicWeight(float v)
{
	float ret = 0.0;
	v = fabs(v);
	float v2 = v*v;
	float v3 = v2*v;
	if (v <= 1.0)
	{
		ret = 1.5 * v3 - 2.5 * v2 + 1.0;
	}
	else if (v > 1.0 && v <= 2.0)
	{
		ret = -0.5 * v3 + 2.5 * v2 - 4.0 * v + 2.0;
	}
	return ret;
}

inline float splineWeight(float v)
{
	float ret = 0.0;
	v = fabs(v);
	float v2 = v*v;
	float v3 = v2*v;
	if (v <= 1.0)
	{
		ret = 0.5 * v3 - v2 + 0.6666666666666667;
	}
	else if (v > 1.0 && v <= 2.0)
	{
		ret = -0.16666666666666667 * v3 + v2 - 2.0 * v + 1.33333333333333333;
	}
	return ret;
}

float CnmiRegistration::getCubicInterpolationValue(float x, float y)
{
	float ret = 0.0;
	int ix = std::floor(x);
	int iy = std::floor(y);
	float wx[4] = { 0 };
	wx[0] = cubicWeight(ix - 1.0 - x);
	wx[1] = cubicWeight(ix - x);
	wx[2] = cubicWeight(ix + 1.0 - x);
	wx[3] = cubicWeight(ix + 2.0 - x);
	float wy[4] = { 0 };
	wy[0] = cubicWeight(iy - 1.0 - y);
	wy[1] = cubicWeight(iy - y);
	wy[2] = cubicWeight(iy + 1.0 - y);
	wy[3] = cubicWeight(iy + 2.0 - y);
	int offy = iy - 1;
	int offx = ix - 1;
	for (int j = 0; j < 4; ++j)
	{
		for (int i = 0; i < 4; ++i)
		{
			ret += m_ptrMoved[(j + offy)*m_width + (i + offx)] * wx[i] * wy[j];
		}
	}
	return ret;
}

float CnmiRegistration::getBsplieInterpolationValue(float x, float y)
{
	float ret = 0.0;
	int ix = std::floor(x);
	int iy = std::floor(y);
	float wx[4] = { 0 };
	wx[0] = splineWeight(ix - 1.0 - x);
	wx[1] = splineWeight(ix - x);
	wx[2] = splineWeight(ix + 1.0 - x);
	wx[3] = splineWeight(ix + 2.0 - x);
	float wy[4] = { 0 };
	wy[0] = splineWeight(iy - 1.0 - y);
	wy[1] = splineWeight(iy - y);
	wy[2] = splineWeight(iy + 1.0 - y);
	wy[3] = splineWeight(iy + 2.0 - y);
	int offy = iy - 1;
	int offx = ix - 1;
	for (int j = 0; j < 4; ++j)
	{
		for (int i = 0; i < 4; ++i)
		{
			ret += m_ptrMoved[(j + offy)*m_width + (i + offx)] * wx[i] * wy[j];
		}
	}
	return ret;
}

std::pair<float, float> CnmiRegistration::getNNInterpolationValueOfGradient(float x, float y)
{
	std::pair<float, float> ret;
	int ix = std::lroundf(x);
	int iy = std::lroundf(y);
	ret = std::make_pair(m_ptrGradientX[iy * m_width + ix], m_ptrGradientY[iy * m_width + ix]);
	return ret;
}

std::pair<float, float> CnmiRegistration::getBilinearInterpolationValueOfGradient(float x, float y)
{
	std::pair<float, float> ret;
	int ix = lround(x);
	int iy = lround(y);
	ret.first = ((1.0 + ix - x) * m_ptrGradientX[iy * m_width + ix] + (x - ix) * m_ptrGradientX[iy * m_width + ix + 1]) * (1.0 + iy - y)
		+ ((1.0 + ix - x) * m_ptrGradientX[(iy + 1) * m_width + ix] + (x - ix) * m_ptrGradientX[(iy + 1) * m_width + ix + 1]) * (y - iy);
	ret.second = ((1.0 + ix - x) * m_ptrGradientY[iy * m_width + ix] + (x - ix) * m_ptrGradientY[iy * m_width + ix + 1]) * (1.0 + iy - y)
		+ ((1.0 + ix - x) * m_ptrGradientY[(iy + 1) * m_width + ix] + (x - ix) * m_ptrGradientY[(iy + 1) * m_width + ix + 1]) * (y - iy);
	return ret;
}

std::pair<float, float> CnmiRegistration::getCubicInterpolationValueOfGradient(float x, float y)
{
	std::pair<float, float> ret;
	int ix = std::floor(x);
	int iy = std::floor(y);
	float wx[4] = { 0 };
	wx[0] = cubicWeight(ix - 1.0 - x);
	wx[1] = cubicWeight(ix - x);
	wx[2] = cubicWeight(ix + 1.0 - x);
	wx[3] = cubicWeight(ix + 2.0 - x);
	float wy[4] = { 0 };
	wy[0] = cubicWeight(iy - 1.0 - y);
	wy[1] = cubicWeight(iy - y);
	wy[2] = cubicWeight(iy + 1.0 - y);
	wy[3] = cubicWeight(iy + 2.0 - y);
	int offy = iy - 1;
	int offx = ix - 1;
	for (int j = 0; j < 4; ++j)
	{
		for (int i = 0; i < 4; ++i)
		{
			ret.first += m_ptrGradientX[(j + offy)*m_width + (i + offx)] * wx[i] * wy[j];
			ret.second += m_ptrGradientY[(j + offy)*m_width + (i + offx)] * wx[i] * wy[j];
		}
	}
	return ret;
}

std::pair<float, float> CnmiRegistration::getBsplieInterpolationValueOfGradient(float x, float y)
{
	std::pair<float, float> ret;
	int ix = std::floor(x);
	int iy = std::floor(y);
	float wx[4] = { 0 };
	wx[0] = splineWeight(ix - 1.0 - x);
	wx[1] = splineWeight(ix - x);
	wx[2] = splineWeight(ix + 1.0 - x);
	wx[3] = splineWeight(ix + 2.0 - x);
	float wy[4] = { 0 };
	wy[0] = splineWeight(iy - 1.0 - y);
	wy[1] = splineWeight(iy - y);
	wy[2] = splineWeight(iy + 1.0 - y);
	wy[3] = splineWeight(iy + 2.0 - y);
	int offy = iy - 1;
	int offx = ix - 1;
	for (int j = 0; j < 4; ++j)
	{
		for (int i = 0; i < 4; ++i)
		{
			ret.first += m_ptrGradientX[(j + offy)*m_width + (i + offx)] * wx[i] * wy[j];
			ret.second += m_ptrGradientY[(j + offy)*m_width + (i + offx)] * wx[i] * wy[j];
		}
	}
	return ret;
}

void CnmiRegistration::getMoveImageminmax(int top, int bottom, int left, int right, float* pmin, float* pmax)
{
	std::function<float(float, float)> f;
	if (m_interpolationType == nn)
	{
		f = std::bind(&CnmiRegistration::getNNInterpolationValue, this, std::placeholders::_1, std::placeholders::_2);
	}
	else if (m_interpolationType == bilinear)
	{
		f = std::bind(&CnmiRegistration::getBilinearInterpolationValue, this, std::placeholders::_1, std::placeholders::_2);
	}
	else if (m_interpolationType == cubic)
	{
		f = std::bind(&CnmiRegistration::getCubicInterpolationValue, this, std::placeholders::_1, std::placeholders::_2);
	}
	else if (m_interpolationType == bsplie)
	{
		f = std::bind(&CnmiRegistration::getBsplieInterpolationValue, this, std::placeholders::_1, std::placeholders::_2);
	}

	if (m_transformType == translation)
	{
		float new_x, new_y;
		for (int j = top; j < bottom; ++j)
		{
			new_y = j + m_transformMatrix[0][1];
			for (int i = left; i < right; ++i)
			{
				new_x = i + m_transformMatrix[0][0];			
				float v = f(new_x, new_y);
				if (v < *pmin) *pmin = v;
				if (v > *pmax) *pmax = v;
			}
		}
	}
	else if (m_transformType == rigid)
	{

	}
	else if (m_transformType == affine)
	{

	}
	else if (m_transformType == perspective)
	{

	}
	else if (m_transformType == freeform)
	{

	}
}

std::vector<double> CnmiRegistration::getMoveImageProbabilityEstimation()
{
	std::vector<double> prob(m_binsCount, 0.0);
	const double p = sqrt(2.0 * R_PI);
	//double h = 0.73;
	double h = 0.4;
	///

	std::function<float(float, float)> fun;
	if (m_interpolationType == nn)
	{
		fun = std::bind(&CnmiRegistration::getNNInterpolationValue, this, std::placeholders::_1, std::placeholders::_2);
	}
	else if (m_interpolationType == bilinear)
	{
		fun = std::bind(&CnmiRegistration::getBilinearInterpolationValue, this, std::placeholders::_1, std::placeholders::_2);
	}
	else if (m_interpolationType == cubic)
	{
		fun = std::bind(&CnmiRegistration::getCubicInterpolationValue, this, std::placeholders::_1, std::placeholders::_2);
	}
	else if (m_interpolationType == bsplie)
	{
		fun = std::bind(&CnmiRegistration::getBsplieInterpolationValue, this, std::placeholders::_1, std::placeholders::_2);
	}
	
	for (int i = 0; i < m_tempSampleCount; ++i)
	{
		float x = m_tempMovedSamplePoint[i].first;
		float y = m_tempMovedSamplePoint[i].second;
		float value = fun(x, y);
		value = value / m_moveBinWidth - m_moveRealMin;
		int val = (int)value;
		if (val < 2.0) val = 2.0;
		if (val > m_binsCount - 3) val = m_binsCount - 3;
		value = val;
		for (int k = 0; k < m_binsCount; ++k)
		{
			//0.73
			//prob[k] += 1.0 / (m_tempSampleCount * p * h) * exp(-0.5 * (k - value)*(k - value) / (h * h));
			//0.4
			prob[k] += 1.0 / (m_tempSampleCount * p * h * 1.084997) * exp(-0.5 * (k - value)*(k - value) / (h * h));
		}
	}

	for (int k = 0; k < m_binsCount; ++k)
	{
		if (prob[k] < 1.0e-50)
		{
			prob[k] = 1.0e-50;
		}
	}

	return prob;
}

std::vector<std::vector<double>> CnmiRegistration::getJointProbabilityEstimation(std::vector<double> &histFix, std::vector<double> &histMove)
{
	std::vector<std::vector<double>> ret;
	
	std::vector<double> tempVec(m_binsCount, 0);
	for (int i = 0; i < m_binsCount; ++i)
	{
		ret.push_back(tempVec);
	}

	histFix = std::vector<double>(m_binsCount, 0);
	histMove = std::vector<double>(m_binsCount, 0);
	std::vector<double> tempHistFix(m_binsCount, 0);
	std::vector<double> tempHistMove(m_binsCount, 0);

	std::function<float(float, float)> fun;
	if (m_interpolationType == nn)
	{
		fun = std::bind(&CnmiRegistration::getNNInterpolationValue, this, std::placeholders::_1, std::placeholders::_2);
	}
	else if (m_interpolationType == bilinear)
	{
		fun = std::bind(&CnmiRegistration::getBilinearInterpolationValue, this, std::placeholders::_1, std::placeholders::_2);
	}
	else if (m_interpolationType == cubic)
	{
		fun = std::bind(&CnmiRegistration::getCubicInterpolationValue, this, std::placeholders::_1, std::placeholders::_2);
	}
	else if (m_interpolationType == bsplie)
	{
		fun = std::bind(&CnmiRegistration::getBsplieInterpolationValue, this, std::placeholders::_1, std::placeholders::_2);
	}

	const double p = 2.0 * R_PI;
	//double h = 0.73;
	double h = 0.4;

	double h2 = h * h;
	double p2h = sqrt(p) * h * 1.084997;

	//double p_f, p_m;
	double gap;

	int sampletoHistLength = m_binsCount * 2 + 2;
	////
	for (int i = 0; i < m_tempSampleCount; ++i)
	{
		int x = m_tempSamplePoint[i].first;
		int y = m_tempSamplePoint[i].second;
		float x_r = m_tempMovedSamplePoint[i].first;
		float y_r = m_tempMovedSamplePoint[i].second;
		float value_left = m_ptrFix[y*m_width + x];
		float value_right = fun(x_r, y_r);
		value_left = value_left / m_fixBinWidth - m_fixRealMin;
		int val_left = (int)value_left;
		if (val_left < 2) val_left = 2;
		if (val_left > m_binsCount - 3) val_left = m_binsCount - 3;
		value_left = val_left;
		value_right = value_right / m_moveBinWidth - m_moveRealMin;
		int val_right = (int)value_right;
		if (val_right < 2) val_right = 2;
		if (val_right > m_binsCount - 3) val_right = m_binsCount - 3;
		value_right = val_right;

		m_ptrSampletoHistValue[i * sampletoHistLength] = value_left;
		m_ptrSampletoHistValue[i * sampletoHistLength + 1] = value_right;
		 
		//for (int k = 0; k < m_binsCount; ++k)
		//{
		//	gap = fabs(k - value_left);
		//	if (gap < 3.0)
		//	{
		//		p_f = 1.0 / p2h * exp(-0.5 * gap * gap / h2);
		//		histFix[k] += p_f;
		//	}
		//	gap = fabs(k - value_right);
		//	if (gap < 3.0)
		//	{
		//		p_m = 1.0 / p2h * exp(-0.5 * gap * gap / h2);
		//		histMove[k] += p_m;
		//	}
		//	for (int m = 0; m < m_binsCount; ++m)
		//	{
		//		//0.73
		//		//ret[k][m] += 1.0 / (m_tempSampleCount * p * h * h) * exp(-0.5 * ((k - value_left)*(k - value_left) + (m - value_right)*(m - value_right)) / (h*h));
		//		//0.4
		//		//ret[k][m] += 1.0 / (m_tempSampleCount * p * h * h * 1.177219) * exp(-0.5 * ((k - value_left)*(k - value_left) + (m - value_right)*(m - value_right)) / (h*h));
		//		//p_m = 1.0 / p2h * exp(-0.5 * (m - value_right)*(m - value_right) / h2);
		//		//histMove[m] += p_m;
		//		gap = fabs(m - value_right);
		//		if (gap < 3.0)
		//		{
		//			p_m = 1.0 / p2h * exp(-0.5 * gap * gap / h2);
		//			ret[k][m] += p_f * p_m;
		//		}
		//	}
		//}
		tempHistFix = std::vector<double>(m_binsCount, 0);
		tempHistMove = std::vector<double>(m_binsCount, 0);
		for (int k = 0; k < m_binsCount; ++k)
		{
			gap = fabs(k - value_left);
			if (gap < 3.0)
			{
				tempHistFix[k] = 1.0 / p2h * exp(-0.5 * gap * gap / h2);
				//tempHistFix[k] = p_f;
				m_ptrSampletoHistValue[i * sampletoHistLength + 2 + k] = tempHistFix[k];
			}
			gap = fabs(k - value_right);
			if (gap < 3.0)
			{
				tempHistMove[k] = 1.0 / p2h * exp(-0.5 * gap * gap / h2);
				//tempHistMove[k] = p_m;
				m_ptrSampletoHistValue[i * sampletoHistLength + 2 + k + m_binsCount] = tempHistMove[k] * (-1.0 * (k - value_right) / h2);
			}
		}
		for (int k = 0; k < m_binsCount; ++k)
		{
			histFix[k] += tempHistFix[k];
			histMove[k] += tempHistMove[k];
			for (int m = 0; m < m_binsCount; ++m)
			{
				ret[k][m] += tempHistFix[k] * tempHistMove[m];
			}
		}
	}

	for (int k = 0; k < m_binsCount; ++k)
	{
		histFix[k] /= m_tempSampleCount;
		histMove[k] /= m_tempSampleCount;
		if (histFix[k] < 1.0e-50)
		{
			histFix[k] = 1.0e-50;
		}
		if (histMove[k] < 1.0e-50)
		{
			histMove[k] = 1.0e-50;
		}
		for (int m = 0; m < m_binsCount; ++m)
		{
			ret[k][m] /= m_tempSampleCount;
			if (ret[k][m] < 1.0e-50)
			{
				ret[k][m] = 1.0e-50;
			}
		}
	}


	/*double sum = 0.0;
	for (int i = 0; i < bins; ++i)
	{
		for (int j = 0; j < bins; ++j)
		{
			sum += ret[i][j];
		}
	}
	printf("JointProb = %.3f\n", sum);*/
	return ret;
}