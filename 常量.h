#ifndef CONSTANT_H
#define CONSTANT_H
const double PI = 3.1415926535898;//Բ����
const double c = 2.99792458e+08;//����
//GPS
const double GM = 3.986005e+14;//��������
const double we = 7.2921151467e-05;//������ת���ٶ�
const double F = -4.442807633e-10;

    //WGS-84���������,GPS
const double a = 6378136.49;//������
const double flattening = 1 / 298.25642;//����
const double delta = 0.0000001;
const double f1 = 1575420000.0;
const double f2 = 1227600000.0;
const double ave_a = 6371008.7714;//����ƽ���뾶
const double hion = 450000.0;//�����߶�
//����
     //BD��Ƶ
const double cf1 = 1561098000.0;     //BD B1Ƶ��Ƶ��(��λHz)
const double cf2 = 1207140000.0;    // BD B2Ƶ��Ƶ��(��λHz)
const double cf3 = 1268520000.0;    // BD B3Ƶ��Ƶ��(��λHz)
const double GM_C = 3.986004418e+14;//����������������
const double we_c = 7.292115e-05;//����������ת���ٶ�

const long _HOUR_IN_MINUTE = 60L;
const long _MINUTE_IN_SECOND = 60L;
const long	_DAY_IN_HOUR = 24L;
const long	_WEEK_IN_DAY = 7L;
const long	_HOUR_IN_SECOND = _HOUR_IN_MINUTE * _MINUTE_IN_SECOND;
const long	_DAY_IN_SECOND = _DAY_IN_HOUR * _HOUR_IN_SECOND;
const long	_WEEK_IN_SECOND = _WEEK_IN_DAY * _DAY_IN_SECOND;


const double H0 = 0;//�ο����ϵĸ߶�,�¶ȣ���ѹ��ʪ��
const double T0 = 20;
const double P0 = 1013.25;
const double RH0 = 0.50;

const double deltjul = 2400000.5;
const int NN = 6;//6�ײ�ֵ
const int O = 2; //��г��������
//harb
/*const double xx = 5084657.6239;
const double yy = 2670325.3492;
const double zz = -2768480.9482;*/
//jfng
/*const double xx = -2279829.0116731;
const double yy = 5004706.47902805;
const double zz = 3219777.4076273;*/
//cas1
/*const double xx = -901776.134925783;
const double yy = 2409383.24758311;
const double zz = -5816748.52877451;*/

/*const double xx = -2279828.96926868;
const double yy = 5004706.49112338;
const double zz = 3219777.42147052;*/
//TWTF
/*const double xx = -2994428.49808605;
const double yy = 4951309.09172896;
const double zz = 2674496.73881957;*/
//cut0
/*const double xx = -2364337.4846811;
const double yy = 4870285.6358861;
const double zz = -3360809.6083186;*/
//cas1
const double xx = -901776.13651914;
const double yy = 2409383.24496966;
const double zz = -5816748.50561365;
//anmg
/*const double xx = -1270826.8700;
const double yy = 6242631.4460;
const double zz = 307792.4399;*/
#endif