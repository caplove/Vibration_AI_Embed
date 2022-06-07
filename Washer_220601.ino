
/*
****************************************************************************************************************************
업데이트 이력

  Version Modified By   Date      Comments                                 NN Version
  ------- -----------  ---------- ---------------------------------------  ------------------------------------
  1.0.0   Junee        17/11/2021 MFCC(x,y,z),STAT(x,y,z)                  myNeuralNetworkFunction_200NN_1104
  1.0.1   Junee        28/12/2021 MFCC(x,y,z,current),STAT(x,y,z,current)
                                  Calibration 추가 (x,y,z,current)
                                  단위는 G(중력), A(암페어)
  1.1.1   Junee        03/01/2022                                          myNeuralNetworkFunction_100NN_220103
  1.1.2   Junee        04/01/2022 Kalman Filter 적용
  1.1.3   Junee        05/01/2022 80% 
  1.1.4   Junee        07/01/2022 92%                                      myNeuralNetworkFunction_200NN_220103(Ver. 1.0)
  1.1.5   Junee        10/01/2022 OLED 버그 수정 예정                                     
  1.1.6   Junee        11/01/2022 재훈련 데이터 (HW위치 수정)                 myNeuralNetworkFunction_200NN_220105(Ver. 1.1)
  1.2.0   Junee        17/01/2022 재훈련 데이터 (2200+)                      myNeuralNetworkFunction_200NN_220117(Ver. 1.2)
  1.2.1   Junee        18/01/2022 재훈련 데이터 (3000+)
  1.3.0   Junee        18/04/2022 전류미터 5A 기준으로 변경                  산학 / Edge Impulse 입력파일 기준  
  1.4.0   Junee        22/04/2022 전면실장기준, CAL 수정 
  1.4.1   Junee        03/05/2022 Edge Impulse Data loader 포멧 추가
  
  2.0.0   Junee        10/05/2022 데이터 전송 타이밍 향상
  2.1.0   Junee        19/05/2022 데이터 전송 타이밍 전체 수정
  2.1.1   Junee        20/05/2022 IMU AHRS 적용
  2.2.0   Junee        26/05/2022 ei-washer-arduino-1.0.14.zip 반영         edeg impulse : Spectral
  2.3.0   Junee        30/05/2022 ei-washer-arduino-1.0.16.zip 반영         edeg impulse : 1D-CNN
  2.3.1   Junee        01/06/2022 Data Forwarder 센서위치 binary 반영
  2.4.0   Junee        03/06/2022 탁도센서 Regression 기능 추가              <washer_turbidity_inferencing.h>
****************************************************************************************************************************

*/

// Include the Arduino library here (something like your_project_inference.h) 
// In the Arduino IDE see **File > Examples > Your project name - Edge Impulse > Static buffer** to get the exact name
#include <washer_turbidity_inferencing.h>
// #include <washer_inferencing.h>
#include <Arduino_LSM9DS1.h>
// #include "Statistic.h"
#include <defs.h>
#include <types.h>
// #include "rtwtypes.h"
#include <stddef.h>
#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include "SimpleKalmanFilter.h"

#include "max_location_index.h"
#include "max_location_index_terminate.h"
#include "rt_nonfinite.h"

 /* IMU */
// https://github.com/armsp/nano-33-ble-gen/tree/master/imu/arduino_madgwick
#include <MadgwickAHRS.h>
Madgwick filter;
unsigned long micros_per_reading, micros_previous;
float accl_scale, gyro_scale;


 /* OLED */
#include <Arduino.h>
#include <U8g2lib.h>

#ifdef U8X8_HAVE_HW_SPI
#include <SPI.h>
#endif
#ifdef U8X8_HAVE_HW_I2C
#include <Wire.h>
#endif

#define SCL_PLOT 0x03
#define SMPL_DATA_COUNT 512     // 데이터 개수
#define SMPL_FREQ 100           // Hz

#define SENSOR_VOLTAGE 3.3



// #define TDS_OFFSET 790      // 전류미터 ADC offset

/* global */

SimpleKalmanFilter simpleKalmanFilter_x(0.01, 0.01, 0.001);   /* 측정오차,   예측오차,  2processing_noise */
SimpleKalmanFilter simpleKalmanFilter_y(0.01, 0.01, 0.001);   /* mea_error, est_error, processing_noise */
SimpleKalmanFilter simpleKalmanFilter_z(0.01, 0.01, 0.001);   /* mea_error, est_error, processing_noise */
SimpleKalmanFilter simpleKalmanFilter_a(0.01, 0.01, 0.001);   /* mea_error, est_error, processing_noise */
SimpleKalmanFilter simpleKalmanFilter_b(0.01, 0.01, 0.001);   /* mea_error, est_error, processing_noise */
SimpleKalmanFilter simpleKalmanFilter_c(0.01, 0.01, 0.001);   /* mea_error, est_error, processing_noise */
SimpleKalmanFilter simpleKalmanFilter_TDS(0.01, 0.01, 0.001);   /* mea_error, est_error, processing_noise */
SimpleKalmanFilter simpleKalmanFilter_TURBIDITY(0.01, 0.01, 0.001);   /* mea_error, est_error, processing_noise */


boolean Kalman_filter_EN = FALSE;
boolean SENSOR_LOCATION_FRONT = TRUE;

U8G2_SSD1306_128X64_NONAME_1_HW_I2C u8g2(U8G2_R0, /* reset=*/ U8X8_PIN_NONE);  /* OLED */

int analogPin = A3, i=0, count = SMPL_DATA_COUNT;
int max_loop_count = 10;
float cal_x=0,cal_y=0,cal_z=0,cal_TDS=0,cal_TURBIDITY=0,sum_x=0,sum_y=0,sum_z=0,sum_TDS=0,sum_TURBIDITY=0;
float cal_x_fr=0,cal_y_fr=0,cal_z_fr=0;
float cal_a,cal_b,cal_c,sum_a,sum_b,sum_c;

const uint16_t samples = SMPL_DATA_COUNT; //This value MUST ALWAYS be a power of 2
const double samplingFrequency = SMPL_FREQ;
double x_max_freq=0, y_max_freq=0, z_max_freq=0, TDS_max_freq=0;
//double vImag_x[samples] = {0,},vImag_y[samples] = {0,},vImag_z[samples] = {0,},vImag_TDS[samples] = {0,};
double label_out[16];
//double label_out_test[16];
//const double vImag[samples] = {0,};

unsigned long timeCheck,timeCheck_mode2[samples],oldtime_Check;
unsigned long delayBetweenSensing = 4200;  // 평균 간격 16800 us 측정됨.
unsigned long delayBetweenSensing_mode2 = 9000;  // 평균 간격 16800 us 측정됨.

float TDS,TURBIDITY;
float signal_TDS[samples],signal_TURBIDITY[samples];
float x,y,z,a,b,c;
float Accroll, Accpitch,AccYaw,Accroll_LF,Accpitch_LF,AccYaw_LF;  // Visualization


int sum_adc_TDS,sum_adc_TURBIDITY;

// 채널별로 데이터 적용 여부 결정하기
int digit;
int select_channel;
int num_digit(int N);








// Edge Impulse 파라미터 //////////////////////////////////////
#define FREQUENCY_HZ        EI_CLASSIFIER_FREQUENCY
#define INTERVAL_MS         (1000 / (FREQUENCY_HZ + 1))

float TDS_OFFSET = 0,TURBIDITY_OFFSET = 812;  // 

static unsigned long last_interval_ms = 0;  // 함수 종료 후에도 시간 저장위해 static으로 할당 unsigned long은 양수정수 범위로 최대 확보

boolean raw_data_check_en= FALSE;


// to classify 1 frame of data you need EI_CLASSIFIER_DSP_INPUT_FRAME_SIZE values
// window size * 채널수 로 정의할 것
// #define EI_CLASSIFIER_DSP_INPUT_FRAME_SIZE 200 * 8     // 실제 NN 에서의 사용은 200 * 7 
float features[EI_CLASSIFIER_DSP_INPUT_FRAME_SIZE];

// keep track of where we are in the feature array
size_t feature_ix = 0;
///////////////////////////////////////////////////////////////////











/////////////////////////////////////////////////////////////////////
// SETUP
/////////////////////////////////////////////////////////////////////
void setup() {

  // micros_per_reading = 1000000/119;
  // micros_previous = micros();

  /* OLED */
   u8g2.begin();
   u8g2.enableUTF8Print();		// enable UTF8 support for the Arduino print() function
  //u8g2.setFlipMode(0);

   u8g2.setFont(u8g2_font_logisoso20_tf);  // font size : 16,18,20,22, ...
   u8g2.setFontDirection(0);
   u8g2.firstPage();
  do {
    u8g2.setCursor(0, 35);
    u8g2.print("Conn. PC");delay(10);
    u8g2.setCursor(0, 63);
    u8g2.print("ver:2.4.0");		// Korean "Hello World" 
  } while ( u8g2.nextPage() );
  delay(1000);
  



  /* PC.Serial */
  delay(1000);
  Serial.begin(115200);
  // Serial.begin(9600);
  
  while (!Serial);              // Serial Port 연결해야 구동됨. --> OLED만 구동시에는 Disable 할 것
  Serial.println("Serial Port Connected ...");
  Serial.println("Ver: 2.4.0  220603");

  Serial.println("--변경점----------------");
  Serial.println("탁도센서 ");
  // Serial.println("gyro 튀는 노이즈 제거");
  // Serial.println("Edge Impulse 입력파일 기준");
    
                                   





  if (!IMU.begin()) {
  Serial.println("Failed to initialize IMU!");
  while (1);
  }

  IMU.setContinuousMode();  // 팝업 오류 제거

  filter.begin(100);  // IMU

  Serial.println("");Serial.println("");Serial.println("");


  /* OLED */
   u8g2.setFont(u8g2_font_logisoso32_tf);  
   u8g2.setFontDirection(0);
   u8g2.firstPage();
  do {
    u8g2.setCursor(20, 40);
    u8g2.print("MENU");delay(10);
    // u8g2.setCursor(0, 40);
    // u8g2.print("CAL 종료");		// Korean "Hello World" 
  } while ( u8g2.nextPage() );
  delay(1000);

}






/////////////////////////////////////////////////////////////////////
// LOOP
//////////////////////////////////////////////////////////////////////

void loop() {
/* 변수 초기화 */
  double window_length = 4.0,samplingFrequency = 100, lowFreq = 5,highFreq = 50;
//  float x,y,z,a,b,c;
  int loop_count = 0, mode = 0, select = 0;
  int softmax_location[10] = {0,};
  int i=0,j=0,k=0, count=SMPL_DATA_COUNT;

  double double_signal_x_tmp[count]= {0,},double_signal_y_tmp[count]= {0,},double_signal_z_tmp[count]= {0,},double_signal_TDS[count]={0,},double_signal_TURBIDITY[count]={0,};
  //double double_signal_a_tmp[count]= {0,},double_signal_b_tmp[count]= {0,},double_signal_c_tmp[count]= {0,};
  float signal_x_tmp[count] = {0,}, signal_y_tmp[count] = {0,}, signal_z_tmp[count] = {0,};
  //float signal_a_tmp[count] = {0,}, signal_b_tmp[count] = {0,}, signal_c_tmp[count] = {0,};
  float signal_TDS[count] = {0,},signal_TURBIDITY[count] = {0,};

//  float x=0, y=0, z=0;
//  int i=0,j=0,k=0, count=SMPL_DATA_COUNT; 
double softmax_value=0;


  //   /* OLED */
  //  u8g2.setFont(u8g2_font_logisoso22_tf);  
  //  u8g2.setFontDirection(0);
  //  u8g2.firstPage();
  // do {
  //   u8g2.setCursor(0, 40);
  //   u8g2.print("-MENU-");delay(10);
  //   // u8g2.setCursor(0, 40);
  //   // u8g2.print("CAL 종료");		// Korean "Hello World" 
  // } while ( u8g2.nextPage() );
  // delay(1000);






/* SW Switch 에 따라 
1) PC 모드 - 트레이닝용 센서 Data만 보냄.
2) MCU모드
*/
  Serial.println("");Serial.println("");Serial.println("");
  Serial.println("** 모드 선택 **************************** ");
  Serial.println("0) 센서위치선택    ");
  Serial.println("1) Calibration    ");
  Serial.println("2) MATLAB RAW 데이터 전송 모드     ");
  Serial.println("3) MATLAB 특성 데이터 전송 / 판정 모드 ");
  Serial.println("4) Processing Visualization 모드 ");
  Serial.println("5) EDGE IMPULSE 데이터 전송 모드");
  Serial.println("6) EDGE IMPULSE 판정 모드[센서개수선택]");


  while(Serial.available() == 0) {
    }
  select = Serial.parseInt();
  Serial.println(select);
  Serial.println("****************************************");
  Serial.println("");Serial.println("");Serial.println("");




  switch (select)
  {
  
  //////////////////////////////////////////////////////////////////
  case 0 :   /* 센서 위치 */
  //////////////////////////////////////////////////////////////////

  /* OLED */
   u8g2.setFont(u8g2_font_logisoso22_tf);  
   u8g2.setFontDirection(0);
   u8g2.firstPage();
  do {
    u8g2.setCursor(0, 40);
    u8g2.print("0) Location");delay(10);
    // u8g2.setCursor(0, 40);
    // u8g2.print("CAL 종료");		// Korean "Hello World" 
  } while ( u8g2.nextPage() );
  delay(1000);


    Serial.println("               1) 센서위치상단     2) 센서위치정면 ");
    while(Serial.available() == 0) {
    }
    select = Serial.parseInt();
    Serial.println(select);
    
    oldtime_Check=0;

    if (select == 2) {
      SENSOR_LOCATION_FRONT = TRUE;
      }
    else
      {
      SENSOR_LOCATION_FRONT = FALSE; 
      }  
      

      
  break;

  //////////////////////////////////////////////////////////////////
  case 5 :   /* EDGE IMPULSE 데이터 전송*/
  //////////////////////////////////////////////////////////////////
    /* OLED */
   u8g2.setFont(u8g2_font_logisoso22_tf);  
   u8g2.setFontDirection(0);
   u8g2.firstPage();
  do {
    u8g2.setCursor(0, 40);
    u8g2.print("5) EDGE");delay(10);
    // u8g2.setCursor(0, 40);
    // u8g2.print("CAL 종료");		// Korean "Hello World" 
  } while ( u8g2.nextPage() );
  delay(1000);


    Serial.println("               1) EDGE - RAW DATA 전송      2) EDGE - RAW DATA 전송(KALMAN Filter 적용) ");
    while(Serial.available() == 0) {
    }
    select = Serial.parseInt();
    Serial.println(select);
    
    oldtime_Check=0;

    if (select == 2) {
      Kalman_filter_EN = TRUE;
      }
    else
      {
      Kalman_filter_EN = FALSE; 
      }  
      
    while(Serial.available() == 0) {
      edge_data_sendinig_mode();
      }    
      
  break;



  //////////////////////////////////////////////////////////////////
  case 6 :   /* EDGE IMPULSE 판정 모드*/
  //////////////////////////////////////////////////////////////////

  Serial.println(" 적용할 센서 채널을 2진표기하시오. 맨앞에 1추가할 것  1XXXXXXX    ");
  Serial.println(" ex) 111100011    --> 8개 채널중 앞의 3개와 뒤의 2개 활용 ");
  // Serial.println(" ex) 1111         -->           앞의 3개  활용 ");
  Serial.println(" ex) 111100000    --> 8개 채널중 앞의 3개  활용(위와 결과 같음) ");

    while(Serial.available() == 0) {
    }
    select_channel = Serial.parseInt();

    digit = num_digit(select_channel);// 자릿수구하기
    digit = digit-1;


    Serial.print(select_channel);Serial.print("--->");Serial.print(digit);Serial.println(" | 채널정보(적용된 채널은 1 )");


    Serial.println("               1) RAW DATA 판정      2) RAW DATA /w KALMAN Filter 판정 ");
    while(Serial.available() == 0) {
    }
    select = Serial.parseInt();
    Serial.println(select);
    
    oldtime_Check=0;

    if (select == 2) {
      Kalman_filter_EN = TRUE;
      }
    else
      {
      Kalman_filter_EN = FALSE; 
      }  


    while(Serial.available() == 0) {
        loop_EDGE_IMPULSE();
    }
  break;





  //////////////////////////////////////////////////////////////////
  case 4 :   /* Visualization 모드 */
  //////////////////////////////////////////////////////////////////
    /* OLED */
   u8g2.setFont(u8g2_font_logisoso22_tf);  
   u8g2.setFontDirection(0);
   u8g2.firstPage();
  do {
    u8g2.setCursor(0, 40);
    u8g2.print("4) VISUAL");delay(10);
    // u8g2.setCursor(0, 40);
    // u8g2.print("CAL 종료");		// Korean "Hello World" 
  } while ( u8g2.nextPage() );
  delay(1000);


    Serial.println("Processing 프로그램과 연동하시오.");

    Serial.println("               1) RAW DATA 시각화      2) KALMAN Filter 적용 시각화 ");
    while(Serial.available() == 0) {
    }
    select = Serial.parseInt();
    Serial.println(select);
    
    if (select == 2) {
      Kalman_filter_EN = TRUE;
      }
    else
      {
      Kalman_filter_EN = FALSE; 
      }  




    while(Serial.available() == 0) {
    Visualization();
    }
  break;


  //////////////////////////////////////////////////////////////////
  case 1 :   /* Cal모드  */
  //////////////////////////////////////////////////////////////////
    /* OLED */
   u8g2.setFont(u8g2_font_logisoso22_tf);  
   u8g2.setFontDirection(0);
   u8g2.firstPage();
  do {
    u8g2.setCursor(0, 40);
    u8g2.print("1) CAL");delay(10);
    // u8g2.setCursor(0, 40);
    // u8g2.print("CAL 종료");		// Korean "Hello World" 
  } while ( u8g2.nextPage() );
  delay(1000);



 Serial.println("               1) 상단 Cal      2) 전면 Cal ");
    while(Serial.available() == 0) {
    }
    select = Serial.parseInt();
    Serial.println(select);
 
    if (select == 1) {
       // 상단
      calibration_top();

      }
    else if (select == 2)
    {
       // 전면
       calibration_front();
    }



  // 전류 offset (ADC) 확인
  sum_adc_TDS=0;
  for (int i = 0; i < count; i++)
    {
     sum_adc_TDS = sum_adc_TDS + analogRead(A1); 
    }

  TDS_OFFSET = sum_adc_TDS / count;

  // 전류 offset (ADC) 확인
  sum_adc_TURBIDITY=0;
  for (int i = 0; i < count; i++)
    {
     sum_adc_TURBIDITY = sum_adc_TURBIDITY + analogRead(A2); 
    }

  TURBIDITY_OFFSET = sum_adc_TURBIDITY / count;

  Serial.print(" TDS값( 0 expected) / 탁도값 ( 812 expected)  : "); Serial.print(TDS_OFFSET);Serial.print(" / ");Serial.println(TURBIDITY_OFFSET);
  
  break;

  //////////////////////////////////////////////////////////////////
  case 2 :   /* PC모드 - raw 데이터 전송 */
  //////////////////////////////////////////////////////////////////  
    /* OLED */
   u8g2.setFont(u8g2_font_logisoso22_tf);  
   u8g2.setFontDirection(0);
   u8g2.firstPage();
  do {
    u8g2.setCursor(0, 40);
    u8g2.print("2) DATA");delay(10);
    // u8g2.setCursor(0, 40);
    // u8g2.print("CAL 종료");		// Korean "Hello World" 
  } while ( u8g2.nextPage() );
  delay(1000);


    Serial.println("               1) RAW DATA 전송      2) RAW DATA 전송(KALMAN Filter 적용) ");
    while(Serial.available() == 0) {
    }
    select = Serial.parseInt();
    Serial.println(select);
    
    oldtime_Check=0;

    if (select == 2) {
      Kalman_filter_EN = TRUE;
      }
    else
      {
      Kalman_filter_EN = FALSE; 
      }  
      
    while(Serial.available() == 0) {
      data_sendinig_mode();
      }    
      
  break;

  //////////////////////////////////////////////////////////////////   
  case 3 :   /* MCU모드 -10번 측정 Routine */
  //////////////////////////////////////////////////////////////////
    /* OLED */
   u8g2.setFont(u8g2_font_logisoso22_tf);  
   u8g2.setFontDirection(0);
   u8g2.firstPage();
  do {
    u8g2.setCursor(0, 40);
    u8g2.print("3) TEST");delay(10);
    // u8g2.setCursor(0, 40);
    // u8g2.print("CAL 종료");		// Korean "Hello World" 
  } while ( u8g2.nextPage() );
  delay(1000);



    Serial.println("               1) 특성 데이터 전송 모드(횟수 입력 필요함.)  2) 판정 모드");
    while(Serial.available() == 0) {
    }
    select = Serial.parseInt();
    Serial.println(select);
    
    if (select == 2) {
        max_loop_count = 10;
        }
    else
        {
        Serial.println("특성 데이터 전송 횟수를 입력하시오.");
        while(Serial.available() == 0) {
        }
        max_loop_count = Serial.parseInt();
        Serial.println(max_loop_count);
        }
    

    for (int m = 0; m < max_loop_count; m++)      
    { 
//  Serial.println("MCU-MODE test-point");      
    loop_count = m + 1;


    /* 센서 데이터 획득 */
    if (select == 2) {
    Serial.println("** 데이터 획득 시작**************************** ");
    }
    for (int i = 0; i < count; i++)
    {
      // IMU.readGyroscope(a, b, c);
       IMU.readAcceleration(x, y, z);
      
       timeCheck_mode2[i] = micros();   // 프로그램 돌리기 시작한 후 지난 밀리 초 숫자 반환
      delayMicroseconds(delayBetweenSensing_mode2);
      signal_x_tmp[i] =  x - cal_x - cal_x_fr;
      signal_y_tmp[i] =  y - cal_y + cal_y_fr;
      signal_z_tmp[i] =  z + cal_z - cal_z_fr;

      signal_TDS[i] = (SENSOR_VOLTAGE / 1024 ) * (analogRead(A1) - TDS_OFFSET) - cal_TDS;      
      signal_TURBIDITY[i] = (SENSOR_VOLTAGE / 1024 ) * (analogRead(A1) - TURBIDITY_OFFSET) + cal_TURBIDITY;  




      // signal_x_tmp[i] = simpleKalmanFilter_x.updateEstimate(x - cal_x);
      // signal_y_tmp[i] = simpleKalmanFilter_y.updateEstimate(y - cal_y);
      // signal_z_tmp[i] = simpleKalmanFilter_z.updateEstimate(z + cal_z);

      // signal_a_tmp[i] =  a;
      // signal_b_tmp[i] =  b;
      // signal_c_tmp[i] =  c;
      

      double_signal_x_tmp[i] = (double) simpleKalmanFilter_x.updateEstimate(x - cal_x - cal_x_fr);
      double_signal_y_tmp[i] = (double) simpleKalmanFilter_y.updateEstimate(y - cal_y + cal_y_fr);
      double_signal_z_tmp[i] = (double) simpleKalmanFilter_z.updateEstimate(z + cal_z - cal_z_fr);

      double_signal_TDS[i] = (double) simpleKalmanFilter_TDS.updateEstimate( (SENSOR_VOLTAGE/1024) * (analogRead(A1) - TDS_OFFSET) - cal_TDS);    // 3.3V/1024 / 0.185 (V/A)   
      double_signal_TURBIDITY[i] = (double) simpleKalmanFilter_TURBIDITY.updateEstimate((SENSOR_VOLTAGE/1024) * (analogRead(A2) - TURBIDITY_OFFSET) + cal_TURBIDITY);    // 3.3V/1024 / 0.185 (V/A)   


      // double_signal_x_tmp[i] = (double) x - cal_x; // Matlab FFT와 통계 산출용
      // double_signal_y_tmp[i] = (double) y - cal_y; // Matlab FFT와 통계 산출용
      // double_signal_z_tmp[i] = (double) z + cal_z; // Matlab FFT와 통계 산출용
      // double_signal_a_tmp[i] = (double) a - cal_a ; // Matlab FFT와 통계 산출용
      // double_signal_b_tmp[i] = (double) b - cal_b ; // Matlab FFT와 통계 산출용
      // double_signal_c_tmp[i] = (double) c - cal_c ; // Matlab FFT와 통계 산출용
      //double_signal_TDS[i] = (double) 0.0322265 * (analogRead(A1) - TDS_OFFSET) - cal_TDS;    // 3.3V/1024 / 0.1 (V/A)  FFT와 통계 산출용 

      // double_x_Stats.add(double_signal_x_tmp[i]);
      // double_y_Stats.add(double_signal_y_tmp[i]);
      // double_z_Stats.add(double_signal_z_tmp[i]);
      // double_TDS_Stats.add(double_signal_TDS[i]);

      //vImag_x[i] = 0;      vImag_y[i] = 0;      vImag_z[i] = 0;      vImag_TDS[i] = 0;        // FFT용
    }
 
 
    // if (select == 2 ){
    // /* 센서 데이터 Serial Port 출력 */
    //   for (int i = 0; i < count; i++)
    //   {  
    //   Serial.print("Sensor_data:");
    //   Serial.print(double_signal_x_tmp[i]);Serial.print(',');
    //   Serial.print(double_signal_y_tmp[i]);Serial.print(',');
    //   Serial.print(double_signal_z_tmp[i]);Serial.print(',');
    //   // Serial.print(signal_a_tmp[i]);Serial.print(',');
    //   // Serial.print(signal_b_tmp[i]);Serial.print(',');
    //   // Serial.print(signal_c_tmp[i]);Serial.print(',');
    //   Serial.print(double_signal_TDS[i]);Serial.print(',');
    //   Serial.print(timeCheck_mode2[i]-timeCheck_mode2[i-1]);
    //   Serial.write(13);Serial.write(10);  // CR/LF
    //   }
    // }
    // delay(3000);




















// // /* 특징 추출 - FFT  
// //    입력 - 가속, 전류 센서 데이터
// //    출력(4*4) - FFT magnitude, maxFreq  --> features_fft[]
// // */

// // Serial.println("Data:");
// //  PrintVector(double_signal_x_tmp, samples, SCL_TIME);
//   FFT.Windowing(double_signal_x_tmp, samples, FFT_WIN_TYP_HAMMING, FFT_FORWARD);	/* Weigh data */
//  // Serial.println("Weighed data:");
//  // PrintVector(double_signal_x_tmp, samples, SCL_TIME);
//   FFT.Compute(double_signal_x_tmp, vImag_x, samples, FFT_FORWARD); /* Compute FFT */  // 앞은 Real 뒤는 Image --> Real만 계산함.
//  // Serial.println("Computed Real values:");
//  // PrintVector(double_signal_x_tmp, samples, SCL_INDEX);
//  // Serial.println("Computed Imaginary values:");
//  // PrintVector(vImag_x, samples, SCL_INDEX);
//   FFT.ComplexToMagnitude(double_signal_x_tmp, vImag_x, samples); /* Compute magnitudes */   // vImg={0,}
//  // Serial.println("Computed magnitudes:");
// //  PrintVector(double_signal_x_tmp, (samples >> 1), SCL_FREQUENCY);
//   double x_max_freq = FFT.MajorPeak(double_signal_x_tmp, samples, samplingFrequency); /* Peak Freq */
// //  Serial.println(x_max_freq, 2);    // 최대크기의 주파수 위치 소숫점2째 자리까지 표기
// Serial.print("FFT_x:");
// for (int i = 0; i < count; i++)
// {
//   Serial.print(double_signal_x_tmp[i]);
//   Serial.print(",");
// }
// Serial.println("");
// Serial.write(13);Serial.write(10);


// // Serial.println("Data:");
// //  PrintVector(double_signal_x_tmp, samples, SCL_TIME);
//   FFT.Windowing(double_signal_y_tmp, samples, FFT_WIN_TYP_HAMMING, FFT_FORWARD);	/* Weigh data */
//  // Serial.println("Weighed data:");
//  // PrintVector(double_signal_x_tmp, samples, SCL_TIME);
//   FFT.Compute(double_signal_y_tmp, vImag_y, samples, FFT_FORWARD); /* Compute FFT */
//  // Serial.println("Computed Real values:");
//  // PrintVector(double_signal_x_tmp, samples, SCL_INDEX);
//  // Serial.println("Computed Imaginary values:");
//  // PrintVector(vImag_x, samples, SCL_INDEX);
//   FFT.ComplexToMagnitude(double_signal_y_tmp, vImag_y, samples); /* Compute magnitudes */
//  // Serial.println("Computed magnitudes:");
// //  PrintVector(double_signal_x_tmp, (samples >> 1), SCL_FREQUENCY);
//   double y_max_freq = FFT.MajorPeak(double_signal_y_tmp, samples, samplingFrequency); /* Peak Freq */
// //  Serial.println(x_max_freq, 2);    // 최대크기의 주파수 위치 소숫점2째 자리까지 표기
// Serial.print("FFT_y:");
// for (int i = 0; i < count; i++)
// {
//   Serial.print(double_signal_y_tmp[i]);
//   Serial.print(",");
// }
// Serial.println("");
// Serial.write(13);Serial.write(10);


// // Serial.println("Data:");
// //  PrintVector(double_signal_x_tmp, samples, SCL_TIME);
//   FFT.Windowing(double_signal_z_tmp, samples, FFT_WIN_TYP_HAMMING, FFT_FORWARD);	/* Weigh data */
//  // Serial.println("Weighed data:");
//  // PrintVector(double_signal_x_tmp, samples, SCL_TIME);
//   FFT.Compute(double_signal_z_tmp, vImag_z, samples, FFT_FORWARD); /* Compute FFT */
//  // Serial.println("Computed Real values:");
//  // PrintVector(double_signal_x_tmp, samples, SCL_INDEX);
//  // Serial.println("Computed Imaginary values:");
//  // PrintVector(vImag_x, samples, SCL_INDEX);
//   FFT.ComplexToMagnitude(double_signal_z_tmp, vImag_z, samples); /* Compute magnitudes */
//  // Serial.println("Computed magnitudes:");
// //  PrintVector(double_signal_x_tmp, (samples >> 1), SCL_FREQUENCY);
//   double z_max_freq = FFT.MajorPeak(double_signal_z_tmp, samples, samplingFrequency); /* Peak Freq */
// //  Serial.println(x_max_freq, 2);    // 최대크기의 주파수 위치 소숫점2째 자리까지 표기
// Serial.print("FFT_z:");
// for (int i = 0; i < count; i++)
// {
//   Serial.print(double_signal_z_tmp[i]);
//   Serial.print(",");
// }
// Serial.println("");
// Serial.write(13);Serial.write(10);


// //  Serial.println("Data:");
// //  PrintVector(double_signal_TDS, samples, SCL_TIME);
//   FFT.Windowing(double_signal_TDS, samples, FFT_WIN_TYP_HAMMING, FFT_FORWARD);	/* Weigh data */
// //  Serial.println("Weighed data:");
// //  PrintVector(double_signal_TDS, samples, SCL_TIME);
//   FFT.Compute(double_signal_TDS, vImag_TDS, samples, FFT_FORWARD); /* Compute FFT */
// //  Serial.println("Computed Real values:");
// //  PrintVector(double_signal_TDS, samples, SCL_INDEX);
//  // Serial.println("Computed Imaginary values:");
// //  PrintVector(vImag_TDS, samples, SCL_INDEX);
//   FFT.ComplexToMagnitude(double_signal_TDS, vImag_TDS, samples); /* Compute magnitudes */
// //  Serial.println("Computed magnitudes:");
// //  PrintVector(double_signal_TDS, (samples >> 1), SCL_FREQUENCY);
//   double TDS_max_freq = FFT.MajorPeak(double_signal_TDS, samples, samplingFrequency);  /* Peak Freq */
// //  Serial.println(TDS_max_freq, 2);        
//   //while(1);       /* Run Once */

// Serial.print("FFT_TDS:");
// for (int i = 0; i < count; i++)
// {
//   Serial.print(double_signal_TDS[i]);
//   Serial.print(",");
// }
// Serial.println("");
// Serial.write(13);Serial.write(10);

// //delay(2000);

// Serial.println("STATUS: FFT Done");







//   features_x_fft[0] = double_signal_x_tmp[1];
//   features_x_fft[1] = double_signal_x_tmp[2];
//   features_x_fft[2] = double_signal_x_tmp[3];
//   features_x_fft[3] = x_max_freq;
//   features_curr_fft[0] = double_signal_TDS[1];
//   features_curr_fft[1] = double_signal_TDS[2];
//   features_curr_fft[2] = double_signal_TDS[3];
//   features_curr_fft[3] = TDS_max_freq;

// delay(4000);

//   Serial.print("Features:");
//   for (i = 0; i < FEATURES_FFT; i++) {
//     Serial.print(features_x_fft[i]);
//     Serial.print(',');
//   }

//   for (i = 0; i < FEATURES_FFT; i++) {
//     Serial.print(features_curr_fft[i]);
//     Serial.print(',');
//   }






/* 
   특징 추출 - MFCC 
   입력 - 가속, 전류 센서 데이터 
   출력(13*4) - MFCC 계수 -->features_mfcc[] 
*/
// myStylus_EXTRACTOR_mfcc_1013(signal_x_tmp, samplingFrequency, window_length, lowFreq,highFreq, MFCCx);
// myStylus_EXTRACTOR_mfcc_1013(signal_y_tmp, samplingFrequency, window_length, lowFreq,highFreq, MFCCy);
// myStylus_EXTRACTOR_mfcc_1013(signal_z_tmp, samplingFrequency, window_length, lowFreq,highFreq, MFCCz);
// myStylus_EXTRACTOR_mfcc_1013(signal_TDS, samplingFrequency,window_length, lowFreq,highFreq, MFCCTDS);                            

// delay(100);

// Serial.print("MFCCx:");
// for (int i = 0; i < FEATURES_MFCC; i++)
// {
// Serial.print(MFCCx[i]);Serial.print(',');
// }
// Serial.println("");
// Serial.print("MFCCy:");
// for ( int i = 0; i < FEATURES_MFCC; i++)
// {
// Serial.print(MFCCy[i]);Serial.print(',');
// }
// Serial.println("");
// Serial.print("MFCCz:");
// for ( int i = 0; i < FEATURES_MFCC; i++)
// {
// Serial.print(MFCCz[i]);Serial.print(',');
// }
// Serial.println("");
// Serial.print("MFCCTDS:");
// for ( i = 0; i < FEATURES_MFCC; i++)
// {
// Serial.print(MFCCTDS[i]);Serial.print(',');      
// }
// Serial.write(13);Serial.write(10); /* CR/LF */

// Serial.println("STATUS: MFCC Done");












// /* 특징 추출 - 통계 
//    입력 - 가속, 전류 센서 데이터
//    출력(3*4) - 범위, 평균, 분산
// */
// Stat_calculation();


// Serial.print("STATx:");
// for ( int i = 0; i < FEATURES_STAT; i++)
// {
// Serial.print(STATx[i]);Serial.print(',');
// }
// Serial.println("");
// Serial.print("STATy:");
// for ( int i = 0; i < FEATURES_STAT; i++)
// {
// Serial.print(STATy[i]);Serial.print(',');
// }
// Serial.println("");
// Serial.print("STATz:");
// for ( int i = 0; i < FEATURES_STAT; i++)
// {
// Serial.print(STATz[i]);Serial.print(',');
// }
// Serial.println("");
// Serial.print("STATTDS:");
// for ( int i = 0; i < FEATURES_STAT; i++)
// {
// Serial.print(STATTDS[i]);Serial.print(',');
// }
// Serial.println("");
//Serial.println("STATUS: STAT Done");



// // 합치기
// for (int i = 0; i < FEATURES_MFCC; i++)
// {
//   input_features[i                ] = (double) MFCCx[i]; 
//   input_features[i+FEATURES_MFCC*1] = (double) MFCCy[i]; 
//   input_features[i+FEATURES_MFCC*2] = (double) MFCCz[i];
//   input_features[i+FEATURES_MFCC*3] = (double) MFCCTDS[i];
   
// }
// for (int i = 0; i < FEATURES_STAT; i++)
// {
//   input_features[i+FEATURES_MFCC*4                ] =  STATx[i]; 
//   input_features[i+FEATURES_MFCC*4+FEATURES_STAT*1] =  STATy[i]; 
//   input_features[i+FEATURES_MFCC*4+FEATURES_STAT*2] =  STATz[i];
//   input_features[i+FEATURES_MFCC*4+FEATURES_STAT*3] =  STATTDS[i]; 
 
// }


// Serial.print("Features:");
// for (int i = 0; i < FEATURES_MFCC*4+FEATURES_STAT*4; i++)    // 1.0.1 전류추가
// //for (int i = 0; i < FEATURES_MFCC*3+FEATURES_STAT*3; i++)
// {
// Serial.print(input_features[i]);Serial.print(",");
// }
// ;Serial.write(13);Serial.write(10); /* CR/LF */
// delay(10);

// double_x_Stats.clear();
// double_y_Stats.clear();
// double_z_Stats.clear();
// double_TDS_Stats.clear();





if (select == 2)
{

//double input_features_test[48] = {49.94,-2.12,0.16,-0.16,-0.70,-0.20,1.58,1.73,0.90,-0.09,-0.88,0.68,-1.35,49.46,-2.98,-0.43,-0.46,-1.79,0.87,-0.29,-0.06,-2.13,-1.15,0.39,1.43,-1.65,51.36,-2.18,-0.82,0.46,0.81,-0.32,-2.51,1.71,-0.54,-0.56,-0.50,0.28,-0.66,0.00,0.01,0.00,0.00,0.01,0.00,0.01,0.01,0.00};
  /* Initialize function 'myNeuralNetworkFunction_1018' input arguments. */
  /* Initialize function input argument 'x1'. */
  /* Call the entry-point 'myNeuralNetworkFunction_1018'. */

  //  myNeuralNetworkFunction_200NN_220117 (input_features, label_out);
//   myNeuralNetworkFunction1029_200layer (input_features_test, label_out_test);
Serial.print("Softmax: ");
  for (int i = 0; i < 16; i++)
  {
    Serial.print(label_out[i]);
    Serial.print(",");
    // if (softmax_value < label_out[i]) {    
    //   softmax_value = label_out[i];
    //   softmax_location[m] = i;
    //}
  }
  Serial.write(13);Serial.write(10); /* CR/LF */


 /* 최대값 index 위치 찾기 */
//softmax_location[m] = (int) max_location_in_array(label_out);    /* label_out : 1 ~ 16 */ 
softmax_location[m] = (int) max_location_index(label_out);
Serial.print("class: ");delay(10);
Serial.print(softmax_location[m]-1);                              /* m : 0~15로 수정함. */
Serial.write(13);Serial.write(10); /* CR/LF */




   /* oled */
   u8g2.setFont(u8g2_font_logisoso32_tf);  
   u8g2.setFontDirection(0);
   u8g2.firstPage();

  do {
    u8g2.setCursor(0, 15);
  switch (loop_count)       /* m+1 */
    {
  case 1:
       u8g2.print("-");      delay(10);break;
  case 2:
       u8g2.print("--");     delay(10);break;
  case 3:
       u8g2.print("---");    delay(10);break;
  case 4:
       u8g2.print("----");   delay(10);break;
  case 5:
       u8g2.print("-----");   delay(10);break;
  case 6:
       u8g2.print("=----");   delay(10);break;
  case 7:
       u8g2.print("==---");   delay(10);break;
  case 8:
       u8g2.print("===--");   delay(10);break;
  case 9:
       u8g2.print("====-");   delay(10);break;
  case 10:
       u8g2.print("=====");   delay(10);break;
    }
    


  u8g2.setCursor(0, 50);
  switch (softmax_location[m])
  {
  case 1:
  u8g2.print("IDLE");     delay(10);  break;
  case 2:
  u8g2.print("Shirt 1");  delay(10);  break;
  case 3:
  u8g2.print("Shirt 2");  delay(10);  break;
  case 4:
  u8g2.print("Shirt 3");  delay(10);  break;
  case 5:
  u8g2.print("Shirt 4");   delay(10); break;
  case 6:
  u8g2.print("Shirt 5");  delay(10);  break;
  case 7:
  u8g2.print("Suit 1");  delay(10);   break;
  case 8:
  u8g2.print("Suit 2");  delay(10);   break;
  case 9:
  u8g2.print("Suit 3");   delay(10);  break;
  case 10:
  u8g2.print("Suit 4");   delay(10);  break;
  case 11:
  u8g2.print("Suit 5");   delay(10);  break;
  case 12:
  u8g2.print("Coat 1");   delay(10);  break;
  case 13:
  u8g2.print("Coat 2");   delay(10);break;
  case 14:
  u8g2.print("Coat 3");   delay(10); break;
  case 15:
  u8g2.print("Coat 4");   delay(10); break;
  case 16:
  u8g2.print("Coat 5");   delay(10); break;
  }

  } while ( u8g2.nextPage());
  delay(10);




  // 마지막 10회에서 최빈값 출력
  if (loop_count == 10)
  {
  int mode_out = 0;
  mode_out = max_frequent_class(softmax_location);  /* 최빈값 산출 함수 호출*/

  u8g2.setFont(u8g2_font_logisoso22_tf);  
  u8g2.firstPage();
  do {
  u8g2.setCursor(0, 40);delay(10);
  u8g2.print("Class:");delay(10);
 
  u8g2.setCursor(80, 40);delay(10);
  u8g2.print(mode_out-1);delay(10);
  } while ( u8g2.nextPage() );
  delay(500);
  }

  }
    }
  break;

}
}




/////////////////////////////////////////////////////////////////////////////////////////////////////////////
/* 10개중에 최빈값 구하기 */
int max_frequent_class(int classes[10])      // 1,2,2,2,2,1,3,2,2,2 --> mode=2
{
//int n = sizeof(classes) / sizeof(int);
int mode,freq,count=1;
for (int i = 0; i < 10; i++)
{
  freq=1;                       // 빈도수는 초기 1 
for (int j = i + 1; j < 10; j++)
		{
			if (classes[i] == classes[j]) //a[i]를 기준으로 a[j]수 들을 비교한다 ex) a[1] == a[2] , a[1] == a[3] ....
				freq += 1;    //a[i] = a[j]가 같을 경우 freq를 1증가 시킨다
		}
		if (freq >= count)    // 새로운 freq 수가 기존에 저장된 count 수 보다 클 경우
		{
			mode = classes[i];
			count = freq;
		}
}
  return mode;
}
/////////////////////////////////////////////////////////////////////////////////////////////////////////////




/////////////////////////////////////////////////////////////////////////////////////////////////////////////
// /* 최빈값 구하기 */
// int max_frequent_class(int *classes)
// {
// int n = sizeof(classes) / sizeof(int);
// int mode,freq,count=1;
// for (int i = 0; i < n; i++)
// {
//   freq=1;
// for (int j = i + 1; j < n; j++)
// 		{
// 			if (classes[i] == classes[j]) //a[i]를 기준으로 a[j]수 들을 비교한다 ex) a[1] == a[2] , a[1] == a[3] ....
// 				freq += 1;    //a[i] = a[j]가 같을 경우 freq를 1증가 시킨다
// 		}
// 		if (freq >= count)    // 새로운 freq 수가 기존에 저장된 count 수 보다 클 경우
// 		{
// 			mode = classes[i];
// 			count = freq;
// 		}
// }
//   return mode;
// }
/////////////////////////////////////////////////////////////////////////////////////////////////////////////



/////////////////////////////////////////////////////////////////////////////////////////////////////////////
// void PrintVector(double *vData, uint16_t bufferSize, uint8_t scaleType)
// {
//   for (uint16_t i = 0; i < bufferSize; i++)
//   {
//     double abscissa;
//     /* Print abscissa value */
//     switch (scaleType)
//     {
//       case SCL_INDEX:
//         abscissa = (i * 1.0);
// 	break;
//       case SCL_TIME:
//         abscissa = ((i * 1.0) / samplingFrequency);
// 	break;
//       case SCL_FREQUENCY:
//         abscissa = ((i * 1.0 * samplingFrequency) / samples);
// 	break;
//     }
//     Serial.print(abscissa, 6);
//     if(scaleType==SCL_FREQUENCY)
//       Serial.print("Hz");
//     Serial.print(" ");
//     Serial.println(vData[i], 4);
//   }
//   Serial.println();
// }
/////////////////////////////////////////////////////////////////////////////////////////////////////////////



/////////////////////////////////////////////////////////////////////////////////////////////////////////////
// void Stat_calculation() 
// {
// double_x_Stats.count();
// double_x_Stats.minimum();
// double_x_Stats.maximum();
// double_x_Stats.average();
// double_x_Stats.unbiased_stdev();

// double_y_Stats.count();
// double_y_Stats.minimum();
// double_y_Stats.maximum();
// double_y_Stats.average();
// double_y_Stats.unbiased_stdev();

// double_z_Stats.count();
// double_z_Stats.minimum();
// double_z_Stats.maximum();
// double_z_Stats.average();
// double_z_Stats.unbiased_stdev();

// double_TDS_Stats.count();
// double_TDS_Stats.minimum();
// double_TDS_Stats.maximum();
// double_TDS_Stats.average();
// double_TDS_Stats.unbiased_stdev();

// STATx[0]=abs(double_x_Stats.maximum()-double_x_Stats.minimum());
// STATx[1]=double_x_Stats.average();
// STATx[2]=double_x_Stats.unbiased_stdev();
// STATy[0]=abs(double_y_Stats.maximum()-double_y_Stats.minimum());
// STATy[1]=double_y_Stats.average();
// STATy[2]=double_y_Stats.unbiased_stdev();
// STATz[0]=abs(double_z_Stats.maximum()-double_z_Stats.minimum());
// STATz[1]=double_z_Stats.average();
// STATz[2]=double_z_Stats.unbiased_stdev();
// STATTDS[0]=abs(double_TDS_Stats.maximum()-double_TDS_Stats.minimum());
// STATTDS[1]=double_TDS_Stats.average();
// STATTDS[2]=double_TDS_Stats.unbiased_stdev();

// delay(10);
// }
/////////////////////////////////////////////////////////////////////////////////////////////////////////////













/////////////////////////////////////////////////////////////////////////////////////////////////////////////
void data_sendinig_mode()
{
   if (millis() > last_interval_ms + INTERVAL_MS) {        // last_interval_ms 에다가 정확히 10ms 만큼 되면 ....
        last_interval_ms = millis();

        // read sensor data in exactly the same way as in the Data Forwarder example
        IMU.readAcceleration(x, y, z);
        IMU.readGyroscope(a, b, c);
        TDS = ( SENSOR_VOLTAGE / 1024 ) * (analogRead(A1) - TDS_OFFSET);     
        TURBIDITY = ( SENSOR_VOLTAGE / 1024) * (analogRead(A2) - TURBIDITY_OFFSET);   
  Serial.print("Sensor_data:");


if (Kalman_filter_EN) {
  // Serial.print(x-cal_x);Serial.print(',');
  // Serial.print(y-cal_y);Serial.print(',');
  // Serial.print(z+cal_z);Serial.print(',');
  // Serial.print(TDS-cal_TDS);Serial.print(',');
  Serial.print( simpleKalmanFilter_x.updateEstimate(x-cal_x-cal_x_fr));Serial.print(',');
  Serial.print( simpleKalmanFilter_y.updateEstimate(y-cal_y+cal_y_fr));Serial.print(',');
  Serial.print( simpleKalmanFilter_z.updateEstimate(z+cal_z-cal_z_fr));Serial.print(',');  
  Serial.print( simpleKalmanFilter_a.updateEstimate(a-cal_a));Serial.print(',');
  Serial.print( simpleKalmanFilter_b.updateEstimate(b-cal_b));Serial.print(',');
  Serial.print( simpleKalmanFilter_c.updateEstimate(c-cal_c));Serial.print(',');   
  Serial.print( simpleKalmanFilter_TDS.updateEstimate(TDS-cal_TDS));Serial.print(',');
  Serial.print( simpleKalmanFilter_TURBIDITY.updateEstimate(TURBIDITY+cal_TURBIDITY));Serial.print(',');
  
  // Serial.print( simpleKalmanFilter_TDS2.updateEstimate(TDS-cal_TDS));Serial.print(',');
  }
else {
  Serial.print(x-cal_x-cal_x_fr);Serial.print(',');
  Serial.print(y-cal_y+cal_y_fr);Serial.print(',');
  Serial.print(z+cal_z-cal_z_fr);Serial.print(',');
  Serial.print(a-cal_a);Serial.print(',');
  Serial.print(b-cal_b);Serial.print(',');
  Serial.print(c-cal_c);Serial.print(',');
  Serial.print(TDS-cal_TDS);Serial.print(',');
  Serial.print(TURBIDITY+cal_TURBIDITY);Serial.print(',');
}
//Serial.print(oldtimeCheck - timeCheck);
// oldtime_Check = timeCheck;
Serial.write(13);Serial.write(10);  // CR/LF
}
}
/////////////////////////////////////////////////////////////////////////////////////////////////////////////



/////////////////////////////////////////////////////////////////////////////////////////////////////////////
void edge_data_sendinig_mode()   // Edge impulse 데이터 획득용
{


    if (millis() > last_interval_ms + INTERVAL_MS) {        // last_interval_ms 에다가 정확히 10ms 만큼 되면 ....
        last_interval_ms = millis();

        // read sensor data in exactly the same way as in the Data Forwarder example
        IMU.readAcceleration(x, y, z);
        IMU.readGyroscope(a, b, c);
        TDS = ( SENSOR_VOLTAGE / 1024 ) * (analogRead(A1) - TDS_OFFSET);     
        TURBIDITY = ( SENSOR_VOLTAGE / 1024) * (analogRead(A2) - TURBIDITY_OFFSET);    


if (Kalman_filter_EN) {
  // Serial.print(x-cal_x);Serial.print(',');
  // Serial.print(y-cal_y);Serial.print(',');
  // Serial.print(z+cal_z);Serial.print(',');
  // Serial.print(TDS-cal_TDS);Serial.print(',');
  Serial.print( simpleKalmanFilter_x.updateEstimate(x-cal_x-cal_x_fr));Serial.print(',');
  Serial.print( simpleKalmanFilter_y.updateEstimate(y-cal_y+cal_y_fr));Serial.print(',');
  Serial.print( simpleKalmanFilter_z.updateEstimate(z+cal_z-cal_z_fr));Serial.print(',');  
  Serial.print( simpleKalmanFilter_a.updateEstimate(a-cal_a));Serial.print(',');
  Serial.print( simpleKalmanFilter_b.updateEstimate(b-cal_b));Serial.print(',');
  Serial.print( simpleKalmanFilter_c.updateEstimate(c-cal_c));Serial.print(',');   
  Serial.print( simpleKalmanFilter_TDS.updateEstimate(TDS-cal_TDS));Serial.print(',');
  Serial.println( simpleKalmanFilter_TURBIDITY.updateEstimate(TURBIDITY+cal_TURBIDITY));
  
  // Serial.print( simpleKalmanFilter_TDS2.updateEstimate(TDS-cal_TDS));Serial.print(',');
  }
else {

        Serial.print(x-cal_x-cal_x_fr);Serial.print(',');
        Serial.print(y-cal_y+cal_y_fr);Serial.print(',');
        Serial.print(z+cal_z-cal_z_fr);Serial.print(',');
        Serial.print(a-cal_a);Serial.print(',');
        Serial.print(b-cal_b);Serial.print(',');
        Serial.print(c-cal_c);Serial.print(',');
        Serial.print(TDS-cal_TDS);Serial.print(',');
        Serial.println(TURBIDITY+cal_TURBIDITY);


        } 
  }
}
// { 
//   IMU.setContinuousMode();
//   // IMU.readGyroscope(a, b, c);     // 중요 Acc 읽고, Gyro 읽어야 튀는 노이즈 없음. !!!
//   IMU.readAcceleration(x, y, z);
//   IMU.readGyroscope(a, b, c); 
//   // delayMicroseconds(10);
//   // TDS = 0.0322265 * (analogRead(A1) - TDS_OFFSET);   // 포트사양 : 3.3V/1024     센서감도 : 0.1 (V/A)
//   TDS = 0.0174 * (analogRead(A1) - TDS_OFFSET);   // 포트사양 : 3.3V/1024     센서감도 : 0.185 (V/A)

//   timeCheck = micros();   // 프로그램 돌리기 시작한 후 지난 밀리 초 숫자 반환

//   delayMicroseconds(delayBetweenSensing);
//   // Serial.print("Sensor_data:");


// if (Kalman_filter_EN) {
//   // Serial.print(x-cal_x);Serial.print(',');
//   // Serial.print(y-cal_y);Serial.print(',');
//   // Serial.print(z+cal_z);Serial.print(',');
//   // Serial.print(TDS-cal_TDS);Serial.print(',');
//   Serial.print( simpleKalmanFilter_x.updateEstimate(x-cal_x-cal_x_fr));Serial.print(',');
//   Serial.print( simpleKalmanFilter_y.updateEstimate(y-cal_y+cal_y_fr));Serial.print(',');
//   Serial.print( simpleKalmanFilter_z.updateEstimate(z+cal_z-cal_z_fr));Serial.print(',');  
//   Serial.print( simpleKalmanFilter_a.updateEstimate(a-cal_a));Serial.print(',');
//   Serial.print( simpleKalmanFilter_b.updateEstimate(b-cal_b));Serial.print(',');
//   Serial.print( simpleKalmanFilter_c.updateEstimate(c-cal_c));Serial.print(',');   
//   Serial.print( simpleKalmanFilter_TDS.updateEstimate(TDS-cal_TDS)); // Serial.print(',');
//   // Serial.print( simpleKalmanFilter_TDS2.updateEstimate(TDS-cal_TDS));Serial.print(',');
//   }
// else {
//   Serial.print(x-cal_x-cal_x_fr);Serial.print(',');
//   Serial.print(y-cal_y+cal_y_fr);Serial.print(',');
//   Serial.print(z+cal_z-cal_z_fr);Serial.print(',');
//   Serial.print(a-cal_a);Serial.print(',');
//   Serial.print(b-cal_b);Serial.print(',');
//   Serial.print(c-cal_c);Serial.print(',');
//   Serial.print(TDS-cal_TDS); // Serial.print(',');
// }
// //Serial.print(oldtimeCheck - timeCheck);
// oldtime_Check = timeCheck;
// Serial.write(13);Serial.write(10);  // CR/LF
// }
/////////////////////////////////////////////////////////////////////////////////////////////////////////////





void loop_EDGE_IMPULSE() 
{

   int select_channel_enable ;
   int temp_number; 

   temp_number = select_channel;  // 초기값
  //  Serial.print(temp_number);


    if (millis() > last_interval_ms + INTERVAL_MS) {
        last_interval_ms = millis();

        // read sensor data in exactly the same way as in the Data Forwarder example
        IMU.readAcceleration(x, y, z);
        IMU.readGyroscope(a, b, c);
        TDS = ( SENSOR_VOLTAGE / 1024 ) * (analogRead(A1) - TDS_OFFSET);     
        TURBIDITY = ( SENSOR_VOLTAGE / 1024) * (analogRead(A2) - TURBIDITY_OFFSET);    



        if (Kalman_filter_EN) {
            temp_number = temp_number - pow(10,digit);
            // Serial.print(temp_number);Serial.print("   ");


            select_channel_enable = temp_number / pow(10, digit-1);   // 몫
            // Serial.print(select_channel_enable);Serial.print("   ");
            if (select_channel_enable == 1 )
            {
              features[feature_ix++] = simpleKalmanFilter_x.updateEstimate(x-cal_x-cal_x_fr);
              temp_number = temp_number - pow(10, digit-1);             
              // Serial.print(" 1번째 ");
            }
            
            select_channel_enable = temp_number / pow(10, digit-2);
            // Serial.print(select_channel_enable);Serial.print("   ");
            if (select_channel_enable == 1)
            {
            features[feature_ix++] = simpleKalmanFilter_y.updateEstimate(y-cal_y+cal_y_fr);
              temp_number = temp_number - pow(10, digit-2);
              // Serial.print(" 2번째 ");
            }

            select_channel_enable = temp_number / pow(10, digit-3);
            // Serial.print(select_channel_enable);Serial.print("   ");        
            if (select_channel_enable == 1)
            {
            features[feature_ix++] = simpleKalmanFilter_z.updateEstimate(z+cal_z-cal_z_fr);
              temp_number = temp_number - pow(10, digit-3);
              // Serial.print(" 3번째 ");
            }

            select_channel_enable = temp_number / pow(10, digit-4);
            // Serial.print(select_channel_enable);Serial.print("   ");        
            if (select_channel_enable == 1)
            {
            features[feature_ix++] = simpleKalmanFilter_a.updateEstimate(a-cal_a);
              temp_number = temp_number - pow(10, digit-4);
              // Serial.print(" 4번째 ");
            }
            
            select_channel_enable = temp_number / pow(10, digit-5);
            // Serial.print(select_channel_enable);Serial.print("   ");
            if (select_channel_enable == 1)
            {
            features[feature_ix++] = simpleKalmanFilter_b.updateEstimate(b-cal_b);
              temp_number = temp_number - pow(10, digit-5);
              // Serial.print(" 5번째 ");
            }

            select_channel_enable = temp_number / pow(10, digit-6);
            // Serial.print(select_channel_enable);Serial.print("   ");
            if (select_channel_enable == 1)
            {
            features[feature_ix++] = simpleKalmanFilter_c.updateEstimate(c-cal_c);
              temp_number = temp_number - pow(10, digit-6);
              // Serial.print(" 6번째 ");
            }

            select_channel_enable = temp_number / pow(10, digit-7);
            // Serial.print(select_channel_enable);Serial.print("   ");
            if (select_channel_enable == 1)
            {
            features[feature_ix++] = simpleKalmanFilter_TDS.updateEstimate(TDS - cal_TDS);      
              temp_number = temp_number - pow(10, digit-7);
              // Serial.print(" 7번째 ");
            }

            select_channel_enable = temp_number ;   // pow(10, digit-8) --> pow(10,0) ? 
            // Serial.print(select_channel_enable);Serial.print("   ");
            if (select_channel_enable == 1)
            {
            features[feature_ix++] = simpleKalmanFilter_TURBIDITY.updateEstimate(TURBIDITY + cal_TURBIDITY);      
              temp_number = temp_number - pow(10, digit-8);
              // Serial.print(" 8번째 ");
            }

        }


      else
        {
            // https://man-25-1.tistory.com/25
            // 자릿수 구하기 --> 최대 자릿수 --> 해당자릿수가 1이면 출력
            // Serial.print("in_loop:");
            temp_number = temp_number - pow(10,digit);
            // Serial.print(temp_number);Serial.print("   ");


            select_channel_enable = temp_number / pow(10, digit-1);   // 몫
            // Serial.print(select_channel_enable);Serial.print("   ");
            if (select_channel_enable == 1 )
            {
              features[feature_ix++] = x-cal_x-cal_x_fr;
              temp_number = temp_number - pow(10, digit-1);             
              // Serial.print(" 1번째 ");
            }
            
            select_channel_enable = temp_number / pow(10, digit-2);
            // Serial.print(select_channel_enable);Serial.print("   ");
            if (select_channel_enable == 1)
            {
            features[feature_ix++] = y-cal_y+cal_y_fr;
              temp_number = temp_number - pow(10, digit-2);
              // Serial.print(" 2번째 ");
            }

            select_channel_enable = temp_number / pow(10, digit-3);
            // Serial.print(select_channel_enable);Serial.print("   ");        
            if (select_channel_enable == 1)
            {
            features[feature_ix++] = z+cal_z-cal_z_fr;
              temp_number = temp_number - pow(10, digit-3);
              // Serial.print(" 3번째 ");
            }

            select_channel_enable = temp_number / pow(10, digit-4);
            // Serial.print(select_channel_enable);Serial.print("   ");        
            if (select_channel_enable == 1)
            {
            features[feature_ix++] = a-cal_a;
              temp_number = temp_number - pow(10, digit-4);
              // Serial.print(" 4번째 ");
            }
            
            select_channel_enable = temp_number / pow(10, digit-5);
            // Serial.print(select_channel_enable);Serial.print("   ");
            if (select_channel_enable == 1)
            {
            features[feature_ix++] = b-cal_b;
              temp_number = temp_number - pow(10, digit-5);
              // Serial.print(" 5번째 ");
            }

            select_channel_enable = temp_number / pow(10, digit-6);
            // Serial.print(select_channel_enable);Serial.print("   ");
            if (select_channel_enable == 1)
            {
            features[feature_ix++] = c-cal_c;
              temp_number = temp_number - pow(10, digit-6);
              // Serial.print(" 6번째 ");
            }

            select_channel_enable = temp_number / pow(10, digit-7);
            // Serial.print(select_channel_enable);Serial.print("   ");
            if (select_channel_enable == 1)
            {
            features[feature_ix++] = TDS - cal_TDS;      
              temp_number = temp_number - pow(10, digit-7);
              // Serial.print(" 7번째 ");
            }

            select_channel_enable = temp_number ;   // pow(10, digit-8) --> pow(10,0) ? 
            // Serial.print(select_channel_enable);Serial.print("   ");
            if (select_channel_enable == 1)
            {
            features[feature_ix++] = TURBIDITY + cal_TURBIDITY;      
              temp_number = temp_number - pow(10, digit-8);
              // Serial.print(" 8번째 ");
            }
         }

        // Serial.println("");

        // Serial.println(last_interval_ms);




        // features buffer full? then classify!
        if (feature_ix == EI_CLASSIFIER_DSP_INPUT_FRAME_SIZE) {
            ei_impulse_result_t result;

            // create signal from features frame
            signal_t signal;
            numpy::signal_from_buffer(features, EI_CLASSIFIER_DSP_INPUT_FRAME_SIZE, &signal);

            // run classifier
            EI_IMPULSE_ERROR res = run_classifier(&signal, &result, false);
            ei_printf("run_classifier returned: %d\n", res);
            if (res != 0) return;

            // print predictions
            ei_printf("Predictions (DSP: %d ms., Classification: %d ms., Anomaly: %d ms.): \n",
                result.timing.dsp, result.timing.classification, result.timing.anomaly);

            // print the predictions
            for (size_t ix = 0; ix < EI_CLASSIFIER_LABEL_COUNT; ix++) {
                ei_printf("%s:\t%.3f\n", result.classification[ix].label, result.classification[ix].value);
            }
        #if EI_CLASSIFIER_HAS_ANOMALY == 1
            ei_printf("anomaly:\t%.3f\n", result.anomaly);
        #endif

            // reset features frame
            feature_ix = 0;
        }
    }
  }




////////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////////

void ei_printf(const char *format, ...) {
    static char print_buf[1024] = { 0 };

    va_list args;
    va_start(args, format);
    int r = vsnprintf(print_buf, sizeof(print_buf), format, args);
    va_end(args);

    if (r > 0) {
        Serial.write(print_buf);
    }
}








































/////////////////////////////////////////////////////////////////////////////////////////////////////////////
void calibration_top()
{
  
  float  Accroll=0, Accpitch = 0,cal_x_roll_pitch=0,cal_y_roll_pitch=0,cal_z_roll_pitch=0;
  
    cal_x=0,cal_y=0,cal_z=0,cal_TDS=0,cal_TURBIDITY=0,sum_x=0,sum_y=0,sum_z=0,sum_TDS=0,sum_TURBIDITY=0;
    cal_x_fr=0,cal_y_fr=0,cal_z_fr=0;
    cal_a=0,cal_b=0,cal_c=0,sum_a=0,sum_b=0,sum_c=0;
    


    Serial.println("** calibration_top ******************* ");
 
     for (i = 0; i < count*2 ; i++)
    {
    IMU.readAcceleration(x, y, z);
    IMU.readGyroscope(a, b, c);

    TDS = ( SENSOR_VOLTAGE / 1024 ) * (analogRead(A1) - TDS_OFFSET);     
    TURBIDITY = ( SENSOR_VOLTAGE / 1024) * (analogRead(A2) - TURBIDITY_OFFSET);   
          
    sum_x = sum_x + x;      sum_y = sum_y + y;      sum_z = sum_z + z;     sum_TDS = sum_TDS + TDS;
    sum_a = sum_a + a;      sum_b = sum_b + b;      sum_c = sum_c + c;     sum_TURBIDITY = sum_TURBIDITY + TDS;
    delay(1);  
    // cal_x = sum_x / (float) count ;     cal_y = sum_y / (float) count;    cal_z = sum_z / (float) count;
    } 
    cal_x = sum_x / (float) (count*2) ;     cal_y = sum_y / (float) (count*2);    cal_z = 1.0 - (sum_z / (float) (count*2));
    cal_a = sum_a / (float) (count*2) ;     cal_b = sum_b / (float) (count*2);    cal_c = sum_c / (float) (count*2);
    cal_TDS = sum_TDS / (float) (count*2);    cal_TURBIDITY = 10 - sum_TURBIDITY / (float) (count*2);

    Serial.println(" calibration_top 적용할 값 ************");
    Serial.print("  Accx:");Serial.print  (cal_x);
    Serial.print("      Accy:");Serial.print  (cal_y);
    Serial.print("      Accz:");Serial.print  (cal_z);
    Serial.print("      gyrx:");Serial.print  (cal_a);
    Serial.print("      gyry:");Serial.print  (cal_b);
    Serial.print("      gyrz:");Serial.print  (cal_c);
    Serial.print("      TDS:") ;Serial.print  (cal_TDS);
    Serial.print("      탁도:");Serial.println(cal_TURBIDITY);
    Serial.println("");Serial.println("");

    delay(10);
    Serial.println(" Reference 값 ************");
    Serial.println(" Accx:0.00      Accy:0.00      Accz:1.00      gyrx:0.00      gyry:0.00       gyrz:0.00      TDS:0(min)     탁도:10(max)");
    Serial.println("");
    Serial.println(" calibration_top 적용한 값 ************");


    IMU.readAcceleration(x, y, z);  
    IMU.readGyroscope(a, b, c);
 
    // TDS = 0.0322265 * (analogRead(A1) - TDS_OFFSET); 
    TDS = SENSOR_VOLTAGE/1024 * (analogRead(A1) - TDS_OFFSET); 
    TURBIDITY = SENSOR_VOLTAGE/1024 * (analogRead(A1) - TURBIDITY_OFFSET); 

    Serial.print(" Accx:");Serial.print  (x - cal_x );
    Serial.print("      Accy:");Serial.print  (y - cal_y);
    Serial.print("      Accz:");Serial.print  (z + cal_z);
    Serial.print("      gyrx:");Serial.print  (a - cal_a);
    Serial.print("      gyry:");Serial.print  (b - cal_b);
    Serial.print("      gyrz:");Serial.print  (c - cal_c);
    Serial.print("      TDS:") ;Serial.print  (TDS - cal_TDS );
    Serial.print("      탁도:");Serial.println(TURBIDITY + cal_TURBIDITY );
    Serial.println("");Serial.println("");Serial.println("");
    
    delay(10);
  
    /* roll & pitch 계산 */
    // cal_x_roll_pitch = (x-cal_x) / 256;
    // cal_y_roll_pitch = (y-cal_y) / 256;
    // cal_z_roll_pitch = (z+cal_z) / 256;
    cal_x_roll_pitch = x / 256;
    cal_y_roll_pitch = y / 256;
    cal_z_roll_pitch = z / 256;

    Accroll = atan(cal_y_roll_pitch / sqrt(pow(cal_x_roll_pitch, 2) + pow(cal_z_roll_pitch, 2))) * 180 / PI;
    Accpitch = atan(-1 * cal_x_roll_pitch / sqrt(pow(cal_x_roll_pitch, 2) + pow(cal_z_roll_pitch, 2))) * 180 / PI;

    Serial.print("roll(degree) : ");Serial.println(Accroll);
    Serial.print("pitch(degree) : ");Serial.println(Accpitch);


   /* OLED */
   u8g2.setFont(u8g2_font_logisoso18_tf);  
  // u8g2.setFontDirection(0);
   u8g2.firstPage();
  do {
    // u8g2.setCursor(0, 15);
    // u8g2.print("CAL 시작");
    u8g2.setCursor(0, 30);
    u8g2.print("CAL.DONE");delay(10);           
    u8g2.setCursor(0, 63);                  
    u8g2.print(Accroll);delay(10);          
    u8g2.setCursor(70, 63);                 
    u8g2.print(Accpitch);delay(10);		 		  
  } while ( u8g2.nextPage() );
  delay(1000);
}
/////////////////////////////////////////////////////////////////////////////////////////////////////////////


/////////////////////////////////////////////////////////////////////////////////////////////////////////////
void calibration_front()
{
  
  float  Accroll=0, Accpitch = 0,cal_x_roll_pitch=0,cal_y_roll_pitch=0,cal_z_roll_pitch=0;
  
    cal_x=0,cal_y=0,cal_z=0,cal_TDS=0,cal_TURBIDITY=0,sum_x=0,sum_y=0,sum_z=0,sum_TDS=0,sum_TURBIDITY=0;
    cal_x_fr=0,cal_y_fr=0,cal_z_fr=0;
    cal_a=0,cal_b=0,cal_c=0,sum_a=0,sum_b=0,sum_c=0;
    


    Serial.println("** calibration_front ******************* ");
 
     for (i = 0; i < count*2 ; i++)
    {
    IMU.readGyroscope(a, b, c);
    IMU.readAcceleration(x, y, z);

    TDS = ( SENSOR_VOLTAGE / 1024 ) * (analogRead(A1) - TDS_OFFSET);     
    TURBIDITY = ( SENSOR_VOLTAGE / 1024) * (analogRead(A2) - TURBIDITY_OFFSET);   
      
    sum_x = sum_x + x;      sum_y = sum_y + y;      sum_z = sum_z + z;     sum_TDS = sum_TDS + TDS;
    sum_a = sum_a + a;      sum_b = sum_b + b;      sum_c = sum_c + c;     sum_TURBIDITY = sum_TURBIDITY + TURBIDITY;
    delay(1);  
    // cal_x = sum_x / (float) count ;     cal_y = sum_y / (float) count;    cal_z = sum_z / (float) count;
    } 
    cal_x_fr = sum_x / (float) (count*2) ;     cal_y_fr = 1.0 - (sum_y / (float) (count*2));    cal_z_fr = sum_z / (float) (count*2);
    cal_a = sum_a / (float) (count*2) ;     cal_b = sum_b / (float) (count*2);    cal_c = sum_c / (float) (count*2);
    cal_TDS = sum_TDS / (float) (count*2);      cal_TURBIDITY = 10 - sum_TURBIDITY / (float) (count*2);

    Serial.println(" calibration_front 적용할 값 ************");
    Serial.print("  Accx:");Serial.print  (cal_x_fr);
    Serial.print("      Accy:");Serial.print  (cal_y_fr);
    Serial.print("      Accz:");Serial.print  (cal_z_fr);
    Serial.print("      gyrx:");Serial.print  (cal_a);
    Serial.print("      gyry:");Serial.print  (cal_b);
    Serial.print("      gyrz:");Serial.print  (cal_c);
    Serial.print("      TDS:") ;Serial.print  (cal_TDS);
    Serial.print("      탁도:");Serial.println(cal_TURBIDITY);
    Serial.println("");Serial.println("");

    delay(10);
    Serial.println(" Reference 값 ************");    
    Serial.println(" Accx:0.00      Accy:1.00      Accz:0.00      gyrx:0.00      gyry:0.00       gyrz:0.00      TDS:0(min)     탁도:10(max)");
    Serial.println("");

    Serial.println(" calibration_front 적용한 값 ************");

    IMU.readGyroscope(a, b, c);
    IMU.readAcceleration(x, y, z);
    // TDS = 0.0322265 * (analogRead(A1) - TDS_OFFSET); 
    TDS = SENSOR_VOLTAGE / 1024 * (analogRead(A1) - TDS_OFFSET);
    TURBIDITY = SENSOR_VOLTAGE / 1024 * (analogRead(A2) - TURBIDITY_OFFSET);


    Serial.print(" Accx:");Serial.print  (x - cal_x_fr );
    Serial.print("      Accy:");Serial.print  (y + cal_y_fr);
    Serial.print("      Accz:");Serial.print  (z - cal_z_fr);
    Serial.print("      gyrx:");Serial.print  (a - cal_a);
    Serial.print("      gyry:");Serial.print  (b - cal_b);
    Serial.print("      gyrz:");Serial.print  (c - cal_c);
    Serial.print("      TDS:") ;Serial.print  (TDS - cal_TDS );
    Serial.print("      탁도:");Serial.println(TURBIDITY + cal_TURBIDITY );
    Serial.println("");Serial.println("");Serial.println("");
    
    delay(10);
  
    /* roll & pitch 계산 */
    // cal_x_roll_pitch = (x-cal_x) / 256;
    // cal_y_roll_pitch = (y-cal_y) / 256;
    // cal_z_roll_pitch = (z+cal_z) / 256;
    // cal_x_roll_pitch = x / 256;
    // cal_y_roll_pitch = y / 256;
    // cal_z_roll_pitch = z / 256;

    // Accroll = atan(cal_y_roll_pitch / sqrt(pow(cal_x_roll_pitch, 2) + pow(cal_z_roll_pitch, 2))) * 180 / PI;
    // Accpitch = atan(-1 * cal_x_roll_pitch / sqrt(pow(cal_x_roll_pitch, 2) + pow(cal_z_roll_pitch, 2))) * 180 / PI;

    // Serial.print("roll(degree) : ");Serial.println(Accroll);
    // Serial.print("pitch(degree) : ");Serial.println(Accpitch);


   /* OLED */
   u8g2.setFont(u8g2_font_logisoso18_tf);  
  // u8g2.setFontDirection(0);
   u8g2.firstPage();
  do {
    // u8g2.setCursor(0, 15);
    // u8g2.print("CAL 시작");
    u8g2.setCursor(0, 30);
    u8g2.print("CAL.DONE");delay(10);           
    u8g2.setCursor(0, 63);                  
    u8g2.print(Accroll);delay(10);          
    u8g2.setCursor(70, 63);                 
    u8g2.print(Accpitch);delay(10);		 		  
  } while ( u8g2.nextPage() );
  delay(1000);
}
/////////////////////////////////////////////////////////////////////////////////////////////////////////////








/////////////////////////////////////////////////////////////////////////////////////////////////////////////
void Visualization()
{

    float x,y,z,a,b,c,mx, my, mz;
    float roll, pitch, heading;


    if (millis() > last_interval_ms + INTERVAL_MS) {        // last_interval_ms 에다가 정확히 10ms 만큼 되면 ....
    last_interval_ms = millis();

        // read sensor data in exactly the same way as in the Data Forwarder example
        IMU.readAcceleration(x, y, z);
        IMU.readGyroscope(a, b, c);
        IMU.readMagneticField(mx, my, mz);        

        TDS = SENSOR_VOLTAGE / 1024 * (analogRead(A1) - TDS_OFFSET);
        TURBIDITY = SENSOR_VOLTAGE / 1024 * (analogRead(A2) - TURBIDITY_OFFSET);



        if (Kalman_filter_EN) {

        filter.update( simpleKalmanFilter_a.updateEstimate(a-cal_a), simpleKalmanFilter_b.updateEstimate(b-cal_b), simpleKalmanFilter_c.updateEstimate(c-cal_c),
        simpleKalmanFilter_x.updateEstimate(x-cal_x-cal_x_fr), simpleKalmanFilter_y.updateEstimate(y-cal_y+cal_y_fr),
        simpleKalmanFilter_z.updateEstimate(z+cal_z-cal_z_fr), mx, my, mz); //for all 3

        // filter.update(a-cal_a, b-cal_b, c-cal_c, x-cal_x-cal_x_fr, y-cal_y+cal_y_fr, z+cal_z-cal_z_fr, mx, my, mz); //for all 3
        // filter.update(a, b, c, x, y, z, mx, my, mz); //for all 3
        // filter.updateIMU(a, b, c, x, y, z);//only for accl, gyro
        }
        else
        {
        filter.update(a-cal_a, b-cal_b, c-cal_c, x-cal_x-cal_x_fr, y-cal_y+cal_y_fr, z+cal_z-cal_z_fr, mx, my, mz); //for all 3
        }




        roll = filter.getRoll();
        pitch = filter.getPitch();
        heading = filter.getYaw();
        Serial.print("Orientation: ");
        Serial.print(heading);
        Serial.print(" ");
        Serial.print(pitch);
        Serial.print(" ");
        Serial.println(roll);

//        micros_previous = micros_previous + micros_per_reading;

// // IMU.readAcceleration(x, y, z);

//     x = x / 256;
//     y = y / 256;
//     z = z / 256;

//     // x = (x - cal_x) / 256;
//     // y = (y - cal_y) / 256;
//     // z = (z + cal_z) / 256;


//     // http://www.nano-i.com/index.php?mid=salse&m=0&listStyle=webzine&document_srl=1272665
//     // 원하는 yaw 값은 아님. 생각해보기 바람..
//     Accroll  = atan( y / sqrt(pow(x, 2) + pow(z, 2))) * 180 / PI;          
//     Accpitch = atan( x / sqrt(pow(y, 2) + pow(z, 2))) * 180 / PI;
//     AccYaw   = atan(sqrt(pow(x, 2) + pow(y, 2))/ z) * 180 / PI;    
//     // Low-pass filter   이전상태 96% 신규상태 4% 반영
//     Accroll_LF = 0.94 * Accroll_LF + 0.06 * Accroll;
//     Accpitch_LF = 0.94 * Accpitch_LF + 0.06 * Accpitch;
//     AccYaw_LF = 0.94 * AccYaw_LF + 0.06 * AccYaw;
//     Serial.print(Accroll_LF);
//     Serial.print("/");
//     Serial.print(Accpitch_LF);
//     Serial.print("/");
//     Serial.println(AccYaw_LF);

 }    
}
/////////////////////////////////////////////////////////////////////////////////////////////////////////////








/////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Processing 코드
/////////////////////////////////////////////////////////////////////////////////////////////////////////////



// import processing.serial.*;
// Serial myPort;

// float yaw = 0.0;
// float pitch = 0.0;
// float roll = 0.0;

// void setup()
// {
//   size(600, 500, P3D);

//   // if you have only ONE serial port active
//   myPort = new Serial(this, Serial.list()[0], 9600); // if you have only ONE serial port active

//   // if you know the serial port name
//   //myPort = new Serial(this, "COM5:", 9600);                    // Windows
//   //myPort = new Serial(this, "/dev/ttyACM0", 9600);             // Linux
//   //myPort = new Serial(this, "/dev/cu.usbmodem1217321", 9600);  // Mac

//   textSize(16); // set text size
//   textMode(SHAPE); // set text mode to shape
// }

// void draw()
// {
//   serialEvent();  // read and parse incoming serial message
//   background(255); // set background to white
//   lights();

//   translate(width/2, height/2); // set position to centre

//   pushMatrix(); // begin object

//   float c1 = cos(radians(roll));
//   float s1 = sin(radians(roll));
//   float c2 = cos(radians(pitch));
//   float s2 = sin(radians(pitch));
//   float c3 = cos(radians(yaw));
//   float s3 = sin(radians(yaw));
//   applyMatrix( c2*c3, s1*s3+c1*c3*s2, c3*s1*s2-c1*s3, 0,
//                -s2, c1*c2, c2*s1, 0,
//                c2*s3, c1*s2*s3-c3*s1, c1*c3+s1*s2*s3, 0,
//                0, 0, 0, 1);

//   drawArduino();

//   popMatrix(); // end of object

//   // Print values to console
//   print(roll);
//   print("\t");
//   print(pitch);
//   print("\t");
//   print(yaw);
//   println();
// }

// void serialEvent()
// {
//   int newLine = 13; // new line character in ASCII
//   String message;
//   do {
//     message = myPort.readStringUntil(newLine); // read from port until new line
//     if (message != null) {
//       String[] list = split(trim(message), " ");
//       if (list.length >= 4 && list[0].equals("Orientation:")) {
//         yaw = float(list[1]); // convert to float yaw
//         pitch = float(list[2]); // convert to float pitch
//         roll = float(list[3]); // convert to float roll
//       }
//     }
//   } while (message != null);
// }

// void drawArduino()
// {
//   /* function contains shape(s) that are rotated with the IMU */
//   stroke(0, 90, 90); // set outline colour to darker teal
//   fill(0, 130, 130); // set fill colour to lighter teal
//   box(300, 10, 200); // draw Arduino board base shape

//   stroke(0); // set outline colour to black
//   fill(80); // set fill colour to dark grey

//   translate(60, -10, 90); // set position to edge of Arduino box
//   box(170, 20, 10); // draw pin header as box

//   translate(-20, 0, -180); // set position to other edge of Arduino box
//   box(210, 20, 10); // draw other pin header as box
// }
/////////////////////////////////////////////////////////////////////////////////////////////////////////////














/////////////////////////////////////////////////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////////////////////////////////////////////////////
int num_digit(int N) {
   int digit = 0;
   while (N != 0) {
      N = N / 10;
      digit++;
   }
   return digit;
}// 숫자의 자리수 반환함수
/////////////////////////////////////////////////////////////////////////////////////////////////////////////

