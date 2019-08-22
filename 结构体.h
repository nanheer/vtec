#ifndef STRUCT_FILE
#define STRUCT_FILE
#include <vector>
#include"����.h"
using namespace std;
const int MAX_OBS_SATE_NUM =40;
const int MAX_OBS_TYPE = 4;
//ʱ��ṹ��
typedef struct utctime{
	int year,month,day,hour,minute;
	double second;
}UTC;
typedef UTC * ptc;

typedef struct gpstime{
	int weeknum;
	long secondnum;
	double secondfrac;
}GPST;
typedef GPST * pgpst;

typedef struct bdtime{
	int week;
	double second;
}BDT;

typedef BDT * pbdt;
typedef struct juliantime{
	long daynum;
	int secondnum;
	double secondfrac;
}JULIANDAY;
typedef JULIANDAY * pjulian;

typedef struct doytime{
	int yearnum;
	int daynum;
	int secondnum;
	double secondfrac;
}DOY;
typedef DOY * pdoy;

//����ṹ��
typedef struct xyz_coordinate{
	double x;
	double y;
	double z;
}XYZ;
typedef XYZ * pxyz;

typedef struct enu_coordinate{
	double northing;
	double easting;
	double upping;
}ENU;
typedef ENU * penu;

typedef struct enupolar_coordinate{
	double range;
	double azimuth;
	double elevation;;
}ENUPOLAR;
typedef ENUPOLAR * penupolar;

typedef struct blh_coordinate{
	double longitude;
	double latitude;
	double height;
	double sun_lon;
}BLH;
typedef BLH * pblh;
//sp3�ļ��ṹ��
typedef struct single_sp3_ephemeris{//һ������һ��ʱ���xyz����
	UTC utime_n;
	double x, y, z;
}s_sp3_ephe;
typedef struct one_sate_ephemris{//һ����������ʱ���xyz����
	string prn;
	vector<s_sp3_ephe> sate_ephem;
}one_sate_ephem;
typedef struct all_sate_ephemeris{//������������ʱ���xyz����
	one_sate_ephem all_ephem[32];//32������
}all_sate_ephem;
typedef all_sate_ephem* psp3;

/*
//��Ӧvtec�����o�ļ���ȡ�ṹ��
typedef struct observe_value_vtec{//һ������һ��ʱ��Ĺ۲�ֵ
	UTC utime_o;
	double p1, p2, l1, l2;
}obs_value;
typedef struct observe_value{//һ����������ʱ��Ĺ۲�ֵ
	string prn;
	vector<obs_value> obs_data;
}obs_one_sate;
typedef struct all_sate_obs{//������������ʱ��Ĺ۲�ֵ
	ohd obsheaddata;
	obs_one_sate all_sate[32];//32������
}obs;
typedef obs* pobs;*/
//o�ļ��ṹ��
//һ������һ��ʱ����ĸ��۲�ֵ
typedef struct observe_value_vtec{//һ������һ��ʱ��Ĺ۲�ֵ
	double p1, p2, l1, l2;
}obs_value;
//һ��ʱ��Ĺ۲�����
typedef struct observe_value{
	UTC utime_o;
	vector<string> prn_list;//һ��ʱ��۲⵽��������������
	vector<obs_value>one_obs_data;//��Ӧ�Ĺ۲���������Ҳ����
}od;
//o�ļ�ͷ
typedef struct obs_head_data{
	XYZ approx_coordinate;
	ENU antena_height;
	string obsdata_type, station;
}ohd;
//o�ļ�
typedef struct obsfile{
	ohd obsheaddata;
	vector<od> obsdata;
}obs;
typedef obs* pobs; 
//vtec����ṹ��
typedef struct one_vtec_result{
	//double sta_lat, sta_lon;//��վγ�ȡ�����
	UTC rtime;//ʱ��
	double coe[(O+1)*(O+1)];//256����г����ϵ��
	double coe1[(O + 1)*(O + 1)];//256����г����ϵ��
	vector<double> residul;
	//double vtec,vtec1, ionex_vtec,sh_vtec;//vtec��С
	//double mf;//ͶӰ����
	//string prn;
	//double p1, p2, p1s, p2s;
}one_result;
typedef struct result_out{
	vector<one_result> allresult;
}result;
typedef result * pvt;
//ionex�ļ��ṹ��
//����dcb
typedef struct satellitedcb{
	string prn;
	double bias, rms;
}satedcb;
//��վdcb
typedef struct stationdcb{
	string station;
	double bias, rms;
}stadcb;
//vtec����
typedef struct tecmaps{
	UTC vtime;
	double tec_values[71][73];
}vtecmap;
typedef struct ionex_file{
	vector<satedcb> sate_dcb;
	vector<stadcb> sta_dcb;
	vector<vtecmap> vtec_map;
	vector<vtecmap> vtec_rms;
}ionex;
typedef  ionex* pio;
//��dcb�ķ������ϵ�����Ҷ���
/*typedef struct dcb_equations{
	double E0, E1, E2, E3,dp;
}dcb;*/
//dcb�ļ��ṹ��
typedef struct dcb_file{
	string prn,station,obs1,obs2;
	double bias, rms;
}dcb_f;
typedef struct all_dcb_file{
	vector<dcb_f> dcb_val;
}dcb;
typedef  dcb* pdcb;
//��г����ϵ��һ���������һ���ṹ��
typedef struct {
	int degree, order;
	double tec;
}one_line_coeff;
//һ��ʱ�����������
typedef struct{
	vector<one_line_coeff>one_coe;
	double pole_lat, pole_lon;
	UTC jul;
}one_coeff;
//һ��ĵ����ݣ�12����13��ʱ�������
typedef struct {
	vector<one_coeff> coeff;
}sh_file;
typedef sh_file * psh;
typedef struct{
	double legend_val[(O + 1)*(O + 1)];
}one_legend;
typedef struct{
	double le[(O+1)*(O+1)];
}one_legendre;
#endif