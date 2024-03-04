// BinaryConvert.cpp : Этот файл содержит функцию "main". Здесь начинается и заканчивается выполнение программы.
//

#include <iostream>
#include <stdio.h>
#include <fstream>

using namespace std;


// заполняется руками:
const int kUzlov = 902;    // количество узлов
const int Ktr = 840; // количество конечных элементов (треугольников)
const int kt1 = 85; // количество узлов с первым краевым
const int Count6 = Ktr * 4; // количество целых чисел в файле nvtr.dat
const int Count4 = Ktr * 4; // количество целых чисел для записи

// пути к файлам 
string path = "C:/TELMA/PRIMER/";



// Чтение и запись с использованием класса stream по 6 чисел
int ReadFirtsConds(string nameIn, string nameOut)
{
	int i;
	int N6 = 0; // число успешно прочитанных чисел
	int N4 = 0; // число успешно записанных чисел

	cout << endl << " Класс stream, по 6 чисел:  ";

	int rec3;

	ifstream binfile_in(path + nameIn, ios::binary); if (!binfile_in.is_open()) return(-1);
	ofstream txtfile_out(nameOut);	if (!txtfile_out.is_open()) { binfile_in.close(); return(-1); }

	txtfile_out << kt1 << endl;

	for (i = 0; i < kt1; i++)
	{
		binfile_in.read((char*)&rec3, sizeof(rec3));
		N6 += binfile_in.gcount() / sizeof(int);

		/* записываем только первые 4 числа записи в текстовый файл*/
		txtfile_out << rec3 << endl;
		if (txtfile_out.fail()) break; // проверка, что запись 4х чисел прошла успешно

		N4 += 1;
	}

	binfile_in.close();
	txtfile_out.close();

	cout << " прочитано=" << N6 << ", записано=" << N4 << endl;

	if (N6 != kt1) return(-2);
	if (N4 != kt1) return(-3);
	return(0);
}

int ReadElems(string nameIn, string nameOut)
{
	int i;
	int N6 = 0; // число успешно прочитанных чисел
	int N4 = 0; // число успешно записанных чисел

	cout << endl << " Класс stream, по 6 чисел:  ";

	struct elem4 { int ver1, ver2, ver3, nmat; } rec;

	ifstream binfile_in(path + nameIn, ios::binary); if (!binfile_in.is_open()) return(-1);
	ofstream txtfile_out(nameOut);	if (!txtfile_out.is_open()) { binfile_in.close(); return(-1); }

	txtfile_out << Ktr * 2 << endl;

	for (i = 0; i < Ktr * 2; i++)
	{
		binfile_in.read((char*)&rec, sizeof(rec));
		N6 += binfile_in.gcount() / sizeof(int);

		/* записываем только первые 4 числа записи в текстовый файл*/
		txtfile_out << rec.ver1 << " " << rec.ver2 << " " << rec.ver3 << " " << rec.nmat << endl;
		if (txtfile_out.fail()) break; // проверка, что запись 4х чисел прошла успешно

		N4 += 4;
	}

	binfile_in.close();
	txtfile_out.close();

	cout << " прочитано=" << N6 << ", записано=" << N4 << endl;

	if (N6 != Ktr * 8) return(-2);
	if (N4 != Ktr * 8) return(-3);
	return(0);
}

int ReadElemsPr(string nameIn, string nameInMat, string nameOut)
{
	int i;
	int N6 = 0; // число успешно прочитанных чисел
	int N4 = 0; // число успешно записанных чисел

	cout << endl << " Класс stream, по 6 чисел:  ";

	struct elem4 { int ver1, ver2, ver3, ver4, ign1, ign2; } rec;
	int mat;

	ifstream binfile_in(path + nameIn, ios::binary); if (!binfile_in.is_open()) return(-1);
	ifstream binfile_in_mat(path + nameInMat, ios::binary); if (!binfile_in_mat.is_open()) { binfile_in.close(); return(-1); }
	ofstream txtfile_out(nameOut);	if (!txtfile_out.is_open()) { binfile_in.close(); binfile_in_mat.close(); return(-1); }

	txtfile_out << Ktr << endl;

	for (i = 0; i < Ktr; i++)
	{
		binfile_in.read((char*)&rec, sizeof(rec));
		binfile_in_mat.read((char*)&mat, sizeof(mat));
		N6 += binfile_in.gcount() / sizeof(int);
		N6 += binfile_in_mat.gcount() / sizeof(int);

		/* записываем только первые 4 числа записи в текстовый файл*/
		txtfile_out << rec.ver1 << " " << rec.ver2 << " " << rec.ver3 << " " << rec.ver4 << " " << mat << endl;
		if (txtfile_out.fail()) break; // проверка, что запись 4х чисел прошла успешно

		N4 += 5;
	}

	binfile_in.close();
	binfile_in_mat.close();
	txtfile_out.close();

	cout << " прочитано=" << N6 << ", записано=" << N4 << endl;

	if (N6 != Ktr * 7) return(-2);
	if (N4 != Ktr * 5) return(-3);
	return(0);
}

int ReadCoords(string nameIn, string nameOut)
{
	int i;
	int N6 = 0; // число успешно прочитанных чисел
	int N4 = 0; // число успешно записанных чисел

	cout << endl << " Класс stream, по 6 чисел:  ";

	struct coord2 { double ver1, ver2; } rec2;

	ifstream binfile_in(path + nameIn, ios::binary); if (!binfile_in.is_open()) return(-1);
	ofstream txtfile_out(nameOut);	if (!txtfile_out.is_open()) { binfile_in.close(); return(-1); }

	txtfile_out << kUzlov << endl;

	for (i = 0; i < kUzlov; i++)
	{
		binfile_in.read((char*)&rec2, sizeof(rec2));
		N6 += binfile_in.gcount() / sizeof(double);

		/* записываем только первые 4 числа записи в текстовый файл*/
		txtfile_out << rec2.ver1 << " " << rec2.ver2 << endl;
		if (txtfile_out.fail()) break; // проверка, что запись 4х чисел прошла успешно

		N4 += 2;
	}

	binfile_in.close();
	txtfile_out.close();

	cout << " прочитано=" << N6 << ", записано=" << N4 << endl;

	if (N6 != kUzlov * 2) return(-2);
	if (N4 != kUzlov * 2) return(-3);
	return(0);
}


int main()
{


	setlocale(LC_ALL, "");

	string elems = "nvtr_treug";
	string coords = "rz.dat";
	string firstConds = "l1.dat";
	cout << "Начало конвертации" << endl;

	cout << "Первые краевые ... ";
	if (ReadFirtsConds(firstConds, "firstConds.txt") == 0) cout << "удачно.\n";
	else cout << "неудачно.\n";

	cout << "Элементы ... ";
	if (ReadElems(elems, "elemsTr.txt") == 0) cout << "удачно.\n";
	else cout << "неудачно.\n";

	cout << "Элементы ... ";
	if (ReadElemsPr("nvtr.dat", "nvkat2d.dat", "elemsPr.txt") == 0) cout << "удачно.\n";
	else cout << "неудачно.\n";

	cout << "Координаты ... ";
	if (ReadCoords(coords, "coords.txt") == 0) cout << "удачно.\n";
	else cout << "неудачно.\n";

	cout << "Конец.\n";
}
