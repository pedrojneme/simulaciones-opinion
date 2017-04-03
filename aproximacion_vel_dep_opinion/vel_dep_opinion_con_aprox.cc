
#include <stdio.h>
#include <string.h>
#include <fstream>
#include <iomanip>
#include <iostream>
#include<cstdlib> // Esta librer�a es para los n�meros aleatorios
#include <cmath>


using namespace std;
const double Pi = 3.14159265358979323846;


int main(int argc,char *argv[]) {

  //generador numeros aleatorio
  srand(time(NULL));

  //variables de la simulacion
  double Trotation = 10;
  double radio_particula = 1;
  double distancia_interaccion = 4*pow(radio_particula,2);
  double percent_inic_AA, percent_inic_A, percent_inic_B, refizq, refder, omega, fi, vel_media,var;
  int N_particulas, tiempo_final, alpha,modo_distribucion_velocidad;
    //condiciones al ejecutar
  if (argc>1) {
    char * pEnd;

    N_particulas=atoi(argv[1]);		// 1024
    percent_inic_AA=strtod(argv[2], &pEnd);	// 0.3
    percent_inic_A=strtod(argv[3], &pEnd);	// 0.2
    percent_inic_B=strtod(argv[4], &pEnd);	// 0.2
    refizq=strtod(argv[5], &pEnd);		        // = t_ini_disp
    refder=strtod(argv[6], &pEnd);		        // >= 10
    omega=strtod(argv[7], &pEnd);
    fi=strtod(argv[8], &pEnd);
    modo_distribucion_velocidad=atoi(argv[9]); // =0 para modo normal // = 1 para distribucion gaussiana //
    alpha=atoi(argv[10]);			// >= 100000
  }



  //defino matrices fi y omega///////////////////////////////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  double Prob_OMEGA[4][4] = {0};
  double Prob_FI[4][4] = {0};
  for (int j=0;j<4;j++) {for (int ii=0;ii<4;ii++) {Prob_OMEGA[j][ii] = omega;}}
  for (int j=0;j<4;j++) {Prob_OMEGA[j][0] = 0;}
  for (int j=0;j<4;j++) {for (int ii=0;ii<4;ii++) {Prob_FI[j][ii] = fi;}}
  for (int j=0;j<4;j++) {Prob_FI[j][3] = 0;}

  //defino vectores de reflexion///////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  double prob_reflexionDER[4] , prob_reflexionIZQ[4];
  prob_reflexionDER[0]=refder, prob_reflexionDER[1]=refder, prob_reflexionDER[2]=refder, prob_reflexionDER[3]=0.;
	prob_reflexionIZQ[0]=0., prob_reflexionIZQ[1]=refizq, prob_reflexionIZQ[2]=refizq, prob_reflexionIZQ[3]=refizq;


  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ///////////////////////////posicion inicial y opiniones iniciales//////////////////////////////////////////////////////////////////


  char filename_1[200];
  //defino archivo salida:
  sprintf(filename_1,"datos.txt");
  ofstream datos(filename_1);
  datos.precision(8);

  datos << "r="<<refder <<" "<<"l="<<refizq<<" "<<"omega="<<omega<<" "<<"fi="<<fi<<" "<<"modo_distribucion_velocidad="<<modo_distribucion_velocidad << '\n';

  char filename_2[200];
  //defino archivo salida:
  sprintf(filename_2,"taotabla.txt");
  ofstream taotabla(filename_2);
  taotabla.precision(8);


  char filename_4[200];
  //defino archivo salida:
  sprintf(filename_4,"tabla.txt");
  ofstream tabla(filename_4);
  tabla.precision(8);


  char filename_6[200];
  //defino archivo salida:
  sprintf(filename_6,"lambda_interaccion_tabla.txt");
  ofstream lambda_interaccion_tabla(filename_6);
  lambda_interaccion_tabla.precision(8);




  tiempo_final = 18000 ;
  int tiempo_fin_interacciones = 15000;
  int tiempo_inicio = 9000;




  double posicion_x[N_particulas] = {0};
  double posicion_y[N_particulas] = {0};
  double distancia_sq[N_particulas] = {0};
  double angulo[N_particulas] = {0};
  double sizeSys = 250;


  double DeltaT = 0.01;
  double p_cambia_OMEGA,p_cambia_FI,p_cambio_REFIZQ,p_cambio_REFDER;
  double tiempo_medicion = 0;
  double AcuForceX[N_particulas] , AcuForceY[N_particulas] , b[N_particulas], ModForce[N_particulas];
  double cutoff_force=0.075*radio_particula/DeltaT;
  double cutoff_force_sq=pow(cutoff_force, 2.0);
  double tiempo_interaccion[4][4] ;
  double numero_interaccion[4][4] ;
  double tao_interaccion[4][4];
  double lambda_interaccion[4][4];
  double tiempo_entre_interacciones[4]={0};
  double tiempo_total_entre_interacciones[4]={0};
  double lambda[4];
  double tao[4]={0};
  double tiempo_total_interaccion[4]={0};
  double numero_total_interacciones[4]={0};
  double pob_media[4]={0};
  double posicionborde_x[N_particulas] = {0};
  double posicionborde_y[N_particulas] = {0};
  double d_x[N_particulas] = {0};
  double d_y[N_particulas] = {0};
  double d_x_borde[N_particulas] = {0};
  double d_y_borde[N_particulas] = {0};
  double distancia_borde_sq[N_particulas] = {0};



  int grilla[102][102][50];
  int posicion_grilla_x[N_particulas]={0};
  int posicion_grilla_y[N_particulas]={0};
  int interactuaba_antes_[N_particulas][N_particulas];
  int len_grilla[102][102];
  int t = 0;
  int interaccion[N_particulas]={0};
  int en_borde[N_particulas] = {0};
  int cambio_anterior[N_particulas]={0};
  int count[4] = {0};
  int count2 = 0,count3 = 0,count4 = 0;
  int l;
  int h;
  int fin_interaccion;







  double lista_velocidades[] = {0.0015,0.002,0.0025,0.0017,0.0021,0.0023,0.0024,0.0027,0.0012,0.0018};


  for (int iteracion = 0; iteracion < 10; iteracion++) {

    vel_media = lista_velocidades[iteracion];

    for (int i = 0; i < N_particulas; i++) {
      l = 1;
      while(l==1){
        l = 0;
        posicion_x[i] = rand()/((double)RAND_MAX+1)*sizeSys;
        posicion_y[i] = rand()/((double)RAND_MAX+1)*sizeSys;
        angulo[i] = rand()/(double)RAND_MAX*2*Pi;
        for (int j = 1; j < i; j++) {
          distancia_sq[i] = pow(posicion_x[i]-posicion_x[j],2.0) + pow(posicion_y[i]-posicion_y[j],2.0);
          if (distancia_sq[i] < distancia_interaccion){
            l=1;
          }
        }
      }
    }

    // asigno las opiniones iniciales
    int opinion[N_particulas]={0};
    double random ;
    for (int i = 0; i < N_particulas; i++) {
      random = rand()/((double)RAND_MAX+1);
      if (random<percent_inic_AA) {
        opinion[i] = 0;
      }
      else{
        if (random< (percent_inic_AA+percent_inic_A)) {
          opinion[i]=1;
        }
        else{
          if (random< (percent_inic_AA+percent_inic_A+percent_inic_B)) {
            opinion[i]=2;
          }
          else{
            opinion[i]=3;
          }
        }
      }
    }


    /////////////////// asigno la distribucion de velocidades (se utiliza modo_distribucion_velocidad para decidir el modo a utilizar)/////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    double velocidad[N_particulas] = {0};
    //velocidad normal constante:     modo = 0
    if (modo_distribucion_velocidad==0) {
      for (int i = 0; i < N_particulas; i++) {
        velocidad[i] = vel_media;
      }
    }

    //distribucion gaussiana:   modo = 1
    var = 20*vel_media/10;
    double velocidadmed = 0;
    if (modo_distribucion_velocidad==1) {
      double u1,u2;
      for (int i = 0; i < N_particulas; i++) {
        u1 = rand()/((double)RAND_MAX+1);
        u2 = rand()/((double)RAND_MAX+1);
        velocidad[i] = pow(-2*log(u1),0.5)*cos(2*Pi*u2);

        velocidad[i] = velocidad[i]*var + vel_media;
        velocidadmed += velocidad[i];
      }
      velocidadmed = velocidadmed/N_particulas;
      std::cout << velocidadmed << '\n';
    }


    // velocidad dependiente de opinion. modo = 2
    if (modo_distribucion_velocidad==2) {
      for (int i = 0; i < N_particulas; i++) {
        velocidad[i] = vel_media*pow(alpha,opinion[i]);
      }
    }


  ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //preparo para el bucle temporal:reseteo todo.

    for (int j=0;j<4;j++) {for (int ii=0;ii<4;ii++) {Prob_OMEGA[j][ii] = omega;}}
    for (int j=0;j<4;j++) {Prob_OMEGA[j][0] = 0;}
    for (int j=0;j<4;j++) {for (int ii=0;ii<4;ii++) {Prob_FI[j][ii] = fi;}}
    for (int j=0;j<4;j++) {Prob_FI[j][3] = 0;}

    //defino vectores de reflexion///////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    prob_reflexionDER[0]=refder, prob_reflexionDER[1]=refder, prob_reflexionDER[2]=refder, prob_reflexionDER[3]=0.;
    prob_reflexionIZQ[0]=0., prob_reflexionIZQ[1]=refizq, prob_reflexionIZQ[2]=refizq, prob_reflexionIZQ[3]=refizq;





    for (int i = 0; i <= 101; i++) {
      for (int j = 0; j <= 101; j++) {
        len_grilla[i][j] = 0;
        for (int l = 0; l < 50; l++) {
          grilla[i][j][l] = -1;
        }
      }
    }
    for (int i = 0; i < N_particulas; i++) {
      cambio_anterior[i]=-1;
      for (int j = 0; j < N_particulas; j++) {
        interactuaba_antes_[i][j]=-1;
      }
    }


    for (int i = 0; i < 4; i++) {
      pob_media[i]=0;
      tao[i]=0;
      lambda[i]=0;
      tiempo_entre_interacciones[i]=0;
      numero_total_interacciones[i] = 0;
      tiempo_total_entre_interacciones[i] = 0;
      for (int j = 0; j < 4; j++) {
        tiempo_interaccion[i][j]=0;
        numero_interaccion[i][j]=0;
        tao_interaccion[i][j]=0;
        lambda_interaccion[i][j];
      }
    }


    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////comienza el bucle temporal:///////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    t = 0;
    tiempo_medicion = 0;
    fin_interaccion = -1;

    while (t<(tiempo_final*100)) {


      if ((t > tiempo_fin_interacciones*100) && (fin_interaccion < 0)){

        fin_interaccion = 1;


        for (int j=0;j<4;j++) {for (int ii=0;ii<4;ii++) {Prob_OMEGA[j][ii] = 0;}}
        for (int j=0;j<4;j++) {for (int ii=0;ii<4;ii++) {Prob_FI[j][ii] = 0;}}

        for (int i = 0; i < 4; i++) {
          prob_reflexionDER[i]=0;
          prob_reflexionIZQ[i]=0;
        }
      }





      //receteo las grillas a su estado nulo:
      for (int i = 0; i <= 101; i++) {
        for (int j = 0; j <= 101; j++) {
          for (int l = 0; l < len_grilla[i][j] ; l++) {
            grilla[i][j][l] = -1;
          }
          len_grilla[i][j] = 0;
        }
      }

        //primero me fijo las particulas  que estan en el borde para poder aplicar la condicion de contorno periodica:
      for (int i = 0; i < N_particulas; i++) {

        posicionborde_y[i] = posicion_y[i];
        posicionborde_x[i] = posicion_x[i];
        en_borde[i] = 0;

        if (posicion_x[i]<2.5){
          posicionborde_x[i]=sizeSys+posicion_x[i];
          en_borde[i]=1;
        }
        if (posicion_x[i]>(sizeSys-2.5)) {
          posicionborde_x[i] = posicion_x[i]-sizeSys;
          en_borde[i]=1;
        }
        if (posicion_y[i]<2.5){
          posicionborde_y[i]=sizeSys+posicion_y[i];
          en_borde[i]=1;
        }
        if (posicion_y[i]>(sizeSys-2.5)) {
          posicionborde_y[i] = posicion_y[i]-sizeSys;
          en_borde[i]=1;
        }

        // defino las fuerzas de las particulas y defino las pociciones en las grillas:
        interaccion[i]=0;
        AcuForceX[i]=0;
        AcuForceY[i]=0;
        posicion_grilla_x[i] = int(posicion_x[i]/2.5+1);
        posicion_grilla_y[i] = int(posicion_y[i]/2.5+1);
        grilla[posicion_grilla_x[i]][posicion_grilla_y[i]][len_grilla[posicion_grilla_x[i]][posicion_grilla_y[i]]] = i;
        len_grilla[posicion_grilla_x[i]][posicion_grilla_y[i]] += 1;

      }


      //defino las grillas de los costados 0 y 101 que simulan a 100 y 1 respectivamente
      //esto se utiliza para las condiciones de contorno y hacer interactuar las particulas en el borde
      //dado que las particulas no pueden ser asignadas a estas grillas solo generan particulas virtuales para medir interaccion en el borde
      for (int i = 1; i < 101; i++) {
        for(int j = 0; j < len_grilla[100][i] ; j++) {grilla[0][i][j] = grilla[100][i][j];}
        for(int j = 0; j < len_grilla[1][i]; j++) {grilla[101][i][j] = grilla[1][i][j];}
        for(int j = 0; j < len_grilla[i][100]; j++) {grilla[i][0][j] = grilla[i][100][j];}
        for(int j = 0; j < len_grilla[i][1]; j++) {grilla[i][101][j] = grilla[1][i][j];}
      }


      /// comienzo bucle para medir si hay interaccion o no, ademas calculo fuerzas y realizo los cambios de opinion
      for (int i = 0; i < N_particulas; i++) {
        for (int l = -1; l <= 1 ; l++) {
        for (int m = -1; m <= 1 ; m++) {
        for (int j = 0; j < len_grilla[int(posicion_grilla_x[i]+l)][int(posicion_grilla_y[i]+m)]; j++) {

          if (grilla[int(posicion_grilla_x[i]+l)][int(posicion_grilla_y[i]+m)][j] != i){
            h = grilla[int(posicion_grilla_x[i]+l)][int(posicion_grilla_y[i]+m)][j];
            distancia_sq[i] = pow(posicion_x[i]-posicion_x[h],2)+pow(posicion_y[i]-posicion_y[h],2);
            d_x[i] = posicion_x[i]-posicion_x[h];
            d_y[i] = posicion_y[i]-posicion_y[h];

          //si la particula esta en el borde veo si la distancia con h es la distancia normal o la dist en borde
          if (en_borde[i]!=0) {
            distancia_borde_sq[i] = pow(posicionborde_x[i]-posicion_x[h],2)+pow(posicionborde_y[i]-posicion_y[h],2);
            if (distancia_borde_sq[i]<distancia_sq[i]) {
              d_x_borde[i] = posicionborde_x[i]-posicion_x[h];
              d_y_borde[i] = posicionborde_y[i]-posicion_y[h];
              distancia_sq[i] = distancia_borde_sq[i];
              d_x[i] = d_x_borde[i];
              d_y[i] = d_y_borde[i];
            }
          }
          if (distancia_sq[i] <= distancia_interaccion) {     //me fijo si interactuan las particulas

            if (fin_interaccion>0){
              tiempo_interaccion[opinion[h]][opinion[i]] += DeltaT;
              if (interactuaba_antes_[i][h]<0) {
                interactuaba_antes_[i][h]=1;
                numero_interaccion[opinion[h]][opinion[i]] += 1;
                numero_total_interacciones[opinion[i]] += 1;
              }
            }

            interaccion[i] = 1;
            b[i] = (abs(velocidad[i]))*pow((2*radio_particula)-(2*DeltaT*abs(velocidad[i])), 2.0);
            ModForce[i]  = b[i]/(distancia_sq[i]*sqrt(distancia_sq[i]));
            AcuForceX[i] += ModForce[i]*d_x[i];
            AcuForceY[i] += ModForce[i]*d_y[i];

            //cambio de opinion por persuacion
            if (cambio_anterior[i] < 0) {
              p_cambia_OMEGA = DeltaT*Prob_OMEGA[opinion[int(interaccion[i]+0.2)]][int(opinion[i]+0.2)];
              p_cambia_FI = DeltaT*Prob_FI[opinion[int(interaccion[i]+0.2)]][int(opinion[i]+0.2)];
              random = rand()/((double)RAND_MAX+1);
              if (random<p_cambia_OMEGA) {
                opinion[i] += -1;
                cambio_anterior[i]=h;
                cambio_anterior[h]=i;//solo un cambio por encuentro
              }
              else{
                if (random<(p_cambia_OMEGA+p_cambia_FI)) {
                  opinion[i] += 1;
                  cambio_anterior[i]=h;
                  cambio_anterior[h]=i;//solo un cambio por encuentro
                }
              }
           }
          }
          if ((cambio_anterior[i]==h) && (distancia_sq[i]>1.2*distancia_interaccion)){
                cambio_anterior[i]=-1;
          }
          if (distancia_sq[i]>1.2*distancia_interaccion) {
            interactuaba_antes_[i][h]=-1;
          }
        }
        }
        }
      }//fin de buscar si hay interaccion

      if (interaccion[i]==0) {  //si no interactua actualizo por reflexion
          if (fin_interaccion>0){tiempo_total_entre_interacciones[opinion[i]] += DeltaT;}
          p_cambio_REFIZQ = DeltaT * prob_reflexionIZQ[opinion[i]];
          p_cambio_REFDER = DeltaT * prob_reflexionDER[opinion[i]];
          random = rand()/((double)RAND_MAX+1);
          if (random < p_cambio_REFDER) {
            opinion[i] += 1;
            count3 += 1;
          }
          else{
            if (random < (p_cambio_REFDER+p_cambio_REFIZQ)) {
              opinion[i] += -1;
              count3 += 1;
            }
          }
        }
      }

      ///////////////////////////////////////////////////////////////////////
      // actualizo posiciones y aplico las condiciones de contorno etc..////
      ///////////////////////////////////////////////////////////////////////

      for (int i = 0; i < N_particulas; i++) {
        if (modo_distribucion_velocidad == 2) {
         velocidad[i]=vel_media*pow(alpha,opinion[i]) ;
        }
        if (interaccion[i]==0) {
          posicion_x[i] += velocidad[i]*DeltaT*cos(angulo[i]);
          posicion_y[i] += velocidad[i]*DeltaT*sin(angulo[i]);
        }
        if (interaccion[i]!=0) {
          if ((pow(AcuForceX[i],2.0)>cutoff_force_sq) || (pow(AcuForceY[i],2.0)>cutoff_force_sq)) {
                  posicion_x[i] += velocidad[i]*cos(angulo[i])*DeltaT+cutoff_force*DeltaT*AcuForceX[i]/(sqrt(pow(AcuForceX[i], 2.0)+pow(AcuForceY[i], 2.0)));
                  posicion_y[i] += velocidad[i]*sin(angulo[i])*DeltaT+cutoff_force*DeltaT*AcuForceY[i]/(sqrt(pow(AcuForceX[i], 2.0)+pow(AcuForceY[i], 2.0)));
          }
          else {
               posicion_x[i] += velocidad[i]*cos(angulo[i])*DeltaT+AcuForceX[i]*DeltaT;
               posicion_y[i] += velocidad[i]*sin(angulo[i])*DeltaT+AcuForceY[i]*DeltaT;
          }
        }

        //aplico las condiciones de contorno:
        if (posicion_x[i]<=0) {
          posicion_x[i] = sizeSys + posicion_x[i];
        }
        else{
          if (posicion_x[i]>=sizeSys) {
            posicion_x[i] = posicion_x[i]-sizeSys;
          }
        }
        if (posicion_y[i]<=0) {
          posicion_y[i] = sizeSys + posicion_y[i];
        }
        else{
          if (posicion_y[i]>=sizeSys) {
            posicion_y[i] = posicion_y[i]-sizeSys;
          }
        }
        //cambio el angulo de la particula:
        if (rand()/(double)RAND_MAX<(DeltaT/Trotation)) {
            angulo[i] += (rand()/(double)RAND_MAX-0.5)*Pi;
            if (angulo[i]>2*Pi) {
              angulo[i]+=-2*Pi;
            }
            if (angulo[i]<0.0) {
              angulo[i]+= 2*Pi;
            }
        }
        //ahora calculo las poblaciones
        if (t>tiempo_inicio*100) {
          count[opinion[i]] += 1;
        }
      }




      //calculo el valor medio de las poblaciones y printeo lo obtenido:

      if ((t%100000)==0) {
          std::cout << "faltan:"<< int((float(tiempo_final)-float(t)/100)/1000) << " pasos para terminar la iteracion numero "<< iteracion << '\n';
      }


      /*if (false){   //cambiar por true si se decea guardar la dinamica tempporal
          if ((t%100)==0) {
           pobl<<" "<< float(t)/100 << " " << count[0]/float(N_particulas) << " "  << count[1]/float(N_particulas) << " " << count[2]/float(N_particulas) << " " << count[3]/float(N_particulas) << '\n';
          }
       } */


      if ((t>tiempo_inicio*100) && (fin_interaccion<0)) {
        if ((t%100)==0) {
          for (int i = 0; i < 4; i++) {
            pob_media[i] += count[i]/float(N_particulas);
          }
          tiempo_medicion += 1;
        }
      }

      if (((t%100)==0)&&(fin_interaccion>0)) {

          for (int i = 0; i < 4; i++) {
            tao[i] = tiempo_total_entre_interacciones[i]/numero_total_interacciones[i];
            for (int j = 0; j < 4; j++) {
               lambda_interaccion[i][j] = tiempo_interaccion[i][j]/numero_interaccion[i][j];
              }
          }
       }

      for (int i = 0;  i < 4; i++) {count[i]=0;}
      t += 1;
    }

    tabla<<" "<< vel_media<<" "<<pob_media[0]/tiempo_medicion<<" "<<pob_media[1]/tiempo_medicion<<" "<<pob_media[2]/tiempo_medicion<<" "<<pob_media[3]/tiempo_medicion<< '\n';
    taotabla<<" "<<vel_media<<" "<<tao[0]<<" "<<tao[1]<<" "<<tao[2]<<" "<<tao[3]<<'\n';



    lambda_interaccion_tabla<<vel_media;
    for (int i = 0; i < 4; i++) {
      for (int j = 0; j < 4; j++) {
        lambda_interaccion_tabla<<" "<<lambda_interaccion[i][j];
      }
    }
    lambda_interaccion_tabla << '\n';
    std::cout << iteracion << '\n';

  }

}
