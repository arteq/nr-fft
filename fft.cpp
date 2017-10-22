/*  PROGRAM: FFT
 *  WERSJA: 0.1
 *  Autor: Artur GRACKI mailto: arteq(at)arteq(dot)org
 *  OSTATNIA MODYFIKACJA: 2007/04/16 (pon) 12:53:28
 *  KEYWORDS: szybka transformata fouriera, fourier, widmo mocy
 */

#include	<iostream>
#include	<fstream>
#include	<math.h>

// Zamiana elementów miejscami
#define SWAP(a,b) tempr=(a);(a)=(b);(b)=tempr

int ile = 0;

using namespace std;

/* ====================================================================== */

void fft(float data[], unsigned long nn, int isign)
/*
Replaces data[1..2*nn] by its discrete Fourier transform, if isign is input as 1; or replaces
data[1..2*nn] by nn times its inverse discrete Fourier transform, if isign is input as -1.
data is a complex array of length nn or, equivalently, a real array of length 2*nn. nn MUST
be an integer power of 2 (this is not checked for!).*/
{
	unsigned long n,mmax,m,j,istep,i;
	double wtemp,wr,wpr,wpi,wi,theta; 
	float tempr,tempi;
	n=nn << 1;
	j=1;
	for (i=1;i<n;i+=2) { 
		if (j > i) {
			SWAP(data[j],data[i]); 
			SWAP(data[j+1],data[i+1]);
		}
	m=n >> 1;
	while (m >= 2 && j > m) {
		j -= m;
		m >>= 1;
		}
		j += m;
	}

	mmax=2;
	while (n > mmax) { 
		istep=mmax << 1;
		theta=isign*(6.28318530717959/mmax); 
		wtemp=sin(0.5*theta);
		wpr = -2.0*wtemp*wtemp;
		wpi=sin(theta);
		wr=1.0;
		wi=0.0;
		for (m=1;m<mmax;m+=2) { 
			for (i=m;i<=n;i+=istep) {
				j=i+mmax; 
				tempr=wr*data[j]-wi*data[j+1];
				tempi=wr*data[j+1]+wi*data[j];
				data[j]=data[i]-tempr;
				data[j+1]=data[i+1]-tempi;
				data[i] += tempr;
				data[i+1] += tempi;
			}
			wr=(wtemp=wr)*wpr-wi*wpi+wr; 
			wi=wi*wpr+wtemp*wpi+wi;
		}
	mmax=istep;
	}
}

/* ====================================================================== */

// Zapisuje dane do pliku
void zapisz(float dane[], char *nazwa, int nn = 1024)
{
	ofstream fp;
	fp.open(nazwa, ios::out);

	cout << "-> Zapisuje FFT do pliku " << nazwa << endl;
	for(int i=0; i<nn; i+=2)
	{
		// Zapisuje kwadraty współczynników
		if(i <= nn/2)
			fp << (1.0*i)/(2.0*nn) << "\t" << sqrt (dane[i]*dane[i] + dane[i+1] * dane[i+1]) << endl;
		else
			fp << (1.0*nn-i)/(-1.0*nn) << "\t"<< sqrt (dane[i]*dane[i] + dane[i+1] * dane[i+1]) << endl;
	}
	fp.close();
	cout << endl;
}

/* ====================================================================== */

// Zapisuje dane oryginalne
void zapisz_org(float dane[], char *nazwa, int nn = 1024)
{
	ofstream fp;
	fp.open(nazwa, ios::out);

	cout << "-> Zapisuje oryginalne dane do pliku " << nazwa << endl;
	for (int i=0; i < nn; i++)
		fp << i << "\t" << dane[i] << endl;

	fp.close();
}

/* ====================================================================== */

void wczytaj_dane(char plik[], float *dat, int nn)
{
	ifstream fp;
	int i;
	float x,y;

	fp.open(plik, ios::in);
	for (i = 0; i < nn; i++)
	{
		fp >> x;
		fp >> y;

		dat[i] = y;
	}

	fp.close();
	cout << "-> Wczytano dane do pamieci..." << endl;
}

/* ====================================================================== */

// Sprawdza ile będzie danych do wczytania
int rozmiar(char plik[])
{
	ifstream fp;
	int i=0;
	float x,y;

	fp.open(plik, ios::in);
	while ( !fp.eof() )
	{
		fp >> x;
		fp >> y;

		i++;
	}
	fp.close();

	i--;
	cout << "-> Plik " << plik << " zawiera: " << i << " punktow danych" << endl;
	return i;
}

/* ====================================================================== */


int main()
{

	int N = 1024;
	float *tab;
	const float dwa_pi = 6.283185307;

	// Żeby random był bardziej randomowy
	srand(time(0));

 	tab = new float[N];

	// y(t) = sin(wt), w= 2pi/T, T=128, N=1024
	for (int i=0; i<N; i++)
	{
		tab[i]=sin(i * dwa_pi/128.0);
	}
	
	zapisz_org(tab, "sin1_org.dat");
	cout << "-> Licze FFT dla sin #1" << endl;
	fft(tab, N/2, 1);
	zapisz(tab, "sin1_fft.dat");

	/* ====================================================================== */
	
	// y(t) = sin(wt), w= 2pi/T, T=120, N=1024
	for (int i=0; i<N; i++)
	{
		tab[i]=sin(i * dwa_pi/120.0);
	}

	zapisz_org(tab, "sin2_org.dat");
	cout << "-> Licze FFT dla sin #2" << endl;
	fft(tab, N/2, 1);
	zapisz(tab, "sin2_fft.dat");
	
	/* ====================================================================== */

	// y(t) = (-2 ln(x1))^0.5 * cos(2*pi*x2)
	// Szum Gaussa
	for (int i=0; i<N; i++)
	{
		// Losujemy sobie dwa randomiki
		float x1 = drand48();
		float x2 = drand48();

		tab[i] = pow( (-2.0 * log(x1)), 0.5) * cos(dwa_pi * x2); 
	}

	zapisz_org(tab, "gauss1_org.dat");
	cout << "-> Licze FFT dla Guass'a #1" << endl;
	fft(tab, N/2, 1);
	zapisz(tab, "gauss1_fft.dat");

	/* ====================================================================== */
	
	// y(t) = (-2 ln(x1))^0.5 * cos(2*pi*x2)
	// Szum Gaussa dla większej ilości punktów
	N = 8192;
	delete tab;
	tab = new float[N];

	for (int i=0; i<N; i++)
	{
		float x1 = drand48();
		float x2 = drand48();

		tab[i] = pow( (-2.0 * log(x1)), 0.5) * cos(dwa_pi * x2); 
	}

	zapisz_org(tab, "gauss2_org.dat", N);
	cout << "-> Licze FFT dla Guass'a #2" << endl;
	fft(tab, N/2, 1);
	zapisz(tab, "gauss2_fft.dat", N);

	/* ====================================================================== */

	N = rozmiar("zad2.txt");
	delete tab;
	tab = new float[N];

	// Danych dla Lorentza mamy w pliku 16501 ale to nie jest potęga dwójki, więc ciężko będzie z tego
	// policzyć FFT. Żeby dostać jakieś sensowne wyniki obcinam przedział do najbliższej potęgi dwójki
	// czyli w tym przypadku do 2^14=16384, dużo nie tracimy a przynajmniej będzie ładny obrazek :-)
	wczytaj_dane("zad2.txt", tab, 16384);
	cout << "-> Licze FFT dla Lorentza" << endl;
	fft(tab, 16384/2, 1);
	zapisz(tab, "lorentz_fft.dat", 16384);

	delete tab;
	return 0;
}
