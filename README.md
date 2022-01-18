# ImageRegistration_NormalisedMutualInformation
代码复现论文《Parzen-Window Based Normalized Mutual Information for Medical Image Registration》
利用归一化互信息对医学图像进行配准
测试代码如下：
CnmiRegistration cr;
ROIRegion roif, roim;
roif.x = roim.x = roif.y = roim.y = 0;
roif.width = roim.width = roif.height = roim.height = 1024;
roif.y = 240;
roif.height = 540;
cr.initRegistration(pFix_d, pMove_d_rot, 1024, 1024, random, rigid, bilinear, roif, roim, sampleCount, 100, 2.0e-2, binsCount);
std::vector<std::vector<float>> matrix;
std::vector<float> temp;
temp.push_back(0.0);
temp.push_back(0.0);
matrix.push_back(temp);
cr.setTransformMatrix(matrix);
cr.setLearnStepandMoment(0.02, 0.5);
cr.runRegistration();
std::vector<std::vector<float>>ret = cr.getTransformMatrix();
