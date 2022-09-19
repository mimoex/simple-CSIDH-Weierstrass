#include "fp.h"

Point l_to_r_binary(const Point& p, const mpz_class& a, const mpz_class& mod, const mpz_class& n)
{
	string bit;
	size_t bit_size;
	Point result, temp_point;
	
	result.x = 0; result.y = 0;
	result.inf = true;
	
	temp_point = p;

	bit = n.get_str(2);
	bit_size = bit.size();

	for (int i = 0, j = bit_size - 1; i < bit_size; i++, j--) {
		if (bit[j] == '1') {
			result = ec_add(result, temp_point, a, mod);	//add
		}
		temp_point = ec_add(temp_point, temp_point, a, mod);	//dobule
	}
	return result;
}

//Veluの公式
//引数a,bを返す
void Velu(const Point& P, mpz_class& a, mpz_class& b, const mpz_class& order, const mpz_class& mod)
{
	Point Q, result;
	mpz_class Qx2, tQ, tQtemp, uQ, uQtemp, wQ, wQtemp;
	mpz_class t = 0, w = 0;
	int i;

	Q.x = 0; Q.y = 0;
	Q.inf = true;

	for (i = 1; i < order; i++) {
		Q = ec_add(P, Q, a, mod);

		//tQ = 3 * Q.x * Q.x + a;
		pow_fp(Q.x, 2, mod, &Qx2);
		mul_fp(3, Qx2, mod, &tQtemp);
		add_fp(tQtemp, a, mod, &tQ);
		
		//uQ = 2 * Q.y * Q.y;
		pow_fp(Q.y, 2, mod, &tQtemp);
		mul_fp(2, tQtemp, mod, &uQ);

		//wQ = uQ + tQ * Q.x;
		mul_fp(tQ, Q.x, mod, &wQtemp);
		add_fp(uQ, wQtemp, mod, &wQ);

		//t += tQ;
		add_fp(t, tQ, mod, &t);
		//w += wQ;
		add_fp(w, wQ, mod, &w);
	}
	//t *= 5;
	mul_fp(t, 5, mod, &t);
	//w *= 7;
	mul_fp(w, 7, mod, &w);
	//a -= t;
	sub_fp(a, t, mod, &a);
	//b -= w;
	sub_fp(b, w, mod, &b);
}

//鍵の生成
void gen_csidh_key()
{
	int i, j;

	/*** 乱数生成 ***/
	random_device rd;
	default_random_engine eng(rd());
	uniform_int_distribution<int> distr(0, 2);  //0から10までの乱数を生成 //論文では-5から5?

	int a[N], b[N]; //Aさんの鍵,Bさんの鍵

	for (i = 0; i < N; i++) {
		a[i] = distr(eng);
		b[i] = distr(eng);


	}
	cout << "Aさんの鍵:" << ends;
	for (i = 0; i < N; i++) {
		cout << a[i] << " " << ends;
	}
	cout << endl;

	cout << "Bさんの鍵:" << ends;
	for (i = 0; i < N; i++) {
		cout << b[i] << " " << ends;
	}
	cout << endl;


	mpz_class mod;
	mod = "5326738796327623094747867617954605554069371494832722337612446642054009560026576537626892113026381253624626941643949444792662881241621373288942880288065659";

	const mpz_class curve_a = 1, curve_b = 0;
	mpz_class bai, add_temp;

	bool check;

	/*** Aさんのstep1 ***/
	mpz_class A_curve_a = curve_a, A_curve_b = curve_b;

	Point PA;

	for (i = 0; i < N; i++) {
		for (j = 0; j < a[i]; j++) {
			check = 0;
			while (check == 0) {
				PA = gen_point_sqrt(A_curve_a, A_curve_b, mod);
				div_fp((mod + 1), primes[i], mod, &bai);
				PA = l_to_r_binary(PA, A_curve_a, mod, bai);
				if (PA.inf == false) {
					check = 1;
				}
			}
			//Pの位数はl[i]
			cout << " " << primes[i] << ends;
			Velu(PA, A_curve_a, A_curve_b, primes[i], mod);
		}
		cout << ", " << ends;
	}
	cout << "\n\nAさんの公開情報:" << A_curve_a << ", " << A_curve_b << endl;
	/*** Aさんのstep1終了 ***/


	/*** Bさんのstep1 ***/
	mpz_class B_curve_a = curve_a, B_curve_b = curve_b;


	Point PB;

	for (i = 0; i < N; i++) {
		for (j = 0; j < b[i]; j++) {
			check = 0;
			while (check == 0) {
				PB = gen_point_sqrt(B_curve_a, B_curve_b, mod);
				div_fp((mod + 1), primes[i], mod, &bai);
				PB = l_to_r_binary(PB, B_curve_a, mod, bai);
				if (PB.inf == 0) {
					check = 1;
				}
			}
			//Pの位数はl[i]
			cout << " " << primes[i] << ends;
			Velu(PB, B_curve_a, B_curve_b, primes[i], mod);
		}
		cout << ", " << ends;
	}
	cout << "\n\nBさんの公開情報:" << B_curve_a << ", " << B_curve_b << endl;
	/*** Bさんのstep1終了 ***/

	/*** Aさんのstep2 ***/
	//Bさんの公開情報を使う
	mpz_class BA_curve_a = B_curve_a, BA_curve_b = B_curve_b;
	Point PBA;
	for (i = 0; i < N; i++) {
		for (j = 0; j < a[i]; j++) {
			check = 0;
			while (check == 0) {
				PBA = gen_point_sqrt(BA_curve_a, BA_curve_b, mod);
				div_fp((mod + 1), primes[i], mod, &bai);
				PBA = l_to_r_binary(PBA, BA_curve_a, mod, bai);
				if (PBA.inf == 0) {
					check = 1;
				}
			}
			//Pの位数はl[i]
			cout << " " << primes[i] << ends;
			Velu(PBA, BA_curve_a, BA_curve_b, primes[i], mod);
		}
		cout << ", " << ends;
	}
	cout << "\n\nAさんの最終データ:" << BA_curve_a << ", " << BA_curve_b << endl;
	/*** Aさんのstep2終了 ***/

	/*** Bさんのstep2 ***/
	//Aさんの公開情報を使う
	mpz_class AB_curve_a = A_curve_a, AB_curve_b = A_curve_b;
	Point PAB;
	for (i = 0; i < N; i++) {
		for (j = 0; j < b[i]; j++) {
			check = 0;
			while (check == 0) {
				PAB = gen_point_sqrt(AB_curve_a, AB_curve_b, mod);
				div_fp((mod + 1), primes[i], mod, &bai);
				PAB = l_to_r_binary(PAB, AB_curve_a, mod, bai);
				if (PAB.inf == 0) {
					check = 1;
				}
			}
			//Pの位数はl[i]
			cout << " " << primes[i] << ends;
			Velu(PAB, AB_curve_a, AB_curve_b, primes[i], mod);
		}
		cout << ", " << ends;
	}
	cout << "\n\nBさんの最終データ:" << AB_curve_a << ", " << AB_curve_b << endl;
	/*** Aさんのstep2終了 ***/

	/*** A,Bの最終データの比較 ***/
	if (BA_curve_a == AB_curve_a && BA_curve_b == AB_curve_b) {
		cout << "CSIDH OK!!!" << endl;
	}
	else {
		cout << "CSIDH error" << endl;
	}
}


int main()
{
	gen_csidh_key();
}