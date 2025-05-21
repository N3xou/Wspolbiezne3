#include <iostream>
#include <fstream>
#include <cstdint>
#include <stdexcept>
#include <cstring>
#include <windows.h>
#include <thread>
#include <vector>
#include <omp.h>
#include <sstream>
#include <cstdlib>
#include <ctime>
#include <cmath>
#include <random>

#pragma pack(push, 1)
struct BMPHeader {
    uint16_t fileType;      // File type, always 4D42h ("BM")
    uint32_t fileSize;      // Size of the file in bytes
    uint16_t reserved1;     // Reserved, always 0
    uint16_t reserved2;     // Reserved, always 0
    uint32_t offsetData;    // Start position of pixel data (bytes from the beginning of the file)
};

struct BMPInfoHeader {
    uint32_t size;          // Size of this header (40 bytes)
    int32_t width;          // Bitmap width in pixels
    int32_t height;         // Bitmap height in pixels
    uint16_t planes;        // Number of color planes (must be 1)
    uint16_t bitCount;      // Number of bits per pixel
    uint32_t compression;   // Compression type (0 = uncompressed)
    uint32_t sizeImage;     // Image size in bytes
    int32_t xPixelsPerMeter;// Horizontal resolution
    int32_t yPixelsPerMeter;// Vertical resolution
    uint32_t colorsUsed;    // Number of colors in the color palette
    uint32_t colorsImportant;// Number of important colors
};
#pragma pack(pop)
struct Point {
    int x, y;
};

struct BMP {
    BMPHeader header;
    BMPInfoHeader infoHeader;
    uint8_t* data; // Dynamically allocated pixel data
};


BMP readBMP(const std::string& filename) {
    BMP bmp;
    std::ifstream file(filename, std::ios::binary);
    if (!file) {
        throw std::runtime_error("Unable to open BMP file.");
    }

    // Read headers
    file.read(reinterpret_cast<char*>(&bmp.header), sizeof(bmp.header));
    file.read(reinterpret_cast<char*>(&bmp.infoHeader), sizeof(bmp.infoHeader));

    if (bmp.header.fileType != 0x4D42 || bmp.infoHeader.bitCount != 24) {
        throw std::runtime_error("Unsupported BMP format. Only 24-bit BMP is supported.");
    }


    int rowSize = ((bmp.infoHeader.bitCount * bmp.infoHeader.width + 31) / 32) * 4; 
    int dataSize = rowSize * abs(bmp.infoHeader.height);

    // Allocate memory for pixel data
    bmp.data = new uint8_t[dataSize];
    std::memset(bmp.data, 0, dataSize); // Initialize to zero

    // Read pixel data
    file.seekg(bmp.header.offsetData, std::ios::beg);
    file.read(reinterpret_cast<char*>(bmp.data), dataSize);

    return bmp;
}

void saveBMP(const std::string& filename, const BMP& bmp) {
    std::ofstream file(filename, std::ios::binary);
    if (!file) {
        throw std::runtime_error("Unable to save BMP file.");
    }

    file.write(reinterpret_cast<const char*>(&bmp.header), sizeof(bmp.header));
    file.write(reinterpret_cast<const char*>(&bmp.infoHeader), sizeof(bmp.infoHeader));

    int rowSize = ((bmp.infoHeader.bitCount * bmp.infoHeader.width + 31) / 32) * 4;
    int dataSize = rowSize * abs(bmp.infoHeader.height);

    file.write(reinterpret_cast<const char*>(bmp.data), dataSize);
}

// Bresenham
void drawLine(BMP& bmp, int x1, int y1, int x2, int y2, int thickness, uint8_t r, uint8_t g, uint8_t b) {
    int dx = abs(x2 - x1), dy = abs(y2 - y1);
    int sx = (x1 < x2) ? 1 : -1;
    int sy = (y1 < y2) ? 1 : -1;
    int err = dx - dy;

    int rowSize = ((bmp.infoHeader.bitCount * bmp.infoHeader.width + 31) / 32) * 4;

    while (true) {
        for (int tx = -thickness / 2; tx <= thickness / 2; ++tx) {
            for (int ty = -thickness / 2; ty <= thickness / 2; ++ty) {
                int px = x1 + tx;
                int py = y1 + ty;
                if (px >= 0 && px < bmp.infoHeader.width && py >= 0 && py < bmp.infoHeader.height) {
                    int index = py * rowSize + px * 3;
                    bmp.data[index] = b;
                    bmp.data[index + 1] = g;
                    bmp.data[index + 2] = r;
                }
            }
        }

        if (x1 == x2 && y1 == y2) break;
        int e2 = 2 * err;
        if (e2 > -dy) { err -= dy; x1 += sx; }
        if (e2 < dx) { err += dx; y1 += sy; }
    }
}
void drawLineOpenMP(BMP& bmp, int x1, int y1, int x2, int y2, int thickness, uint8_t r, uint8_t g, uint8_t b, int numThreads) {
    std::vector<std::pair<int, int>> points;

    int dx = abs(x2 - x1), dy = abs(y2 - y1);
    int sx = x1 < x2 ? 1 : -1;
    int sy = y1 < y2 ? 1 : -1;
    int err = dx - dy;

    int px = x1, py = y1;
    while (true) {
        points.emplace_back(px, py);
        if (px == x2 && py == y2) break;
        int e2 = 2 * err;
        if (e2 > -dy) { err -= dy; px += sx; }
        if (e2 < dx) { err += dx; py += sy; }
    }

    int rowSize = ((bmp.infoHeader.bitCount * bmp.infoHeader.width + 31) / 32) * 4;


    int chunkSize = points.size() / numThreads;
    int remainder = points.size() % numThreads;

#pragma omp parallel for
    for (int t = 0; t < numThreads; ++t) {
        int start = t * chunkSize;
        int end = (t == numThreads - 1) ? points.size() : start + chunkSize + (t < remainder ? 1 : 0);

        for (int i = start; i < end; ++i) {
            int x = points[i].first;
            int y = points[i].second;
            for (int tx = -thickness / 2; tx <= thickness / 2; ++tx) {
                for (int ty = -thickness / 2; ty <= thickness / 2; ++ty) {
                    int px = x + tx;
                    int py = y + ty;
                    if (px >= 0 && px < bmp.infoHeader.width && py >= 0 && py < bmp.infoHeader.height) {
                        int index = py * rowSize + px * 3;
                        bmp.data[index] = b;
                        bmp.data[index + 1] = g;
                        bmp.data[index + 2] = r;
                    }
                }
            }
        }
    }
}

void drawLineThreaded(BMP& bmp, int x1, int y1, int x2, int y2, int thickness, uint8_t r, uint8_t g, uint8_t b, int numThreads) {
    std::vector<std::pair<int, int>> points;

    int dx = abs(x2 - x1), dy = abs(y2 - y1);
    int sx = x1 < x2 ? 1 : -1;
    int sy = y1 < y2 ? 1 : -1;
    int err = dx - dy;

    int px = x1, py = y1;
    while (true) {
        points.emplace_back(px, py);
        if (px == x2 && py == y2) break;
        int e2 = 2 * err;
        if (e2 > -dy) { err -= dy; px += sx; }
        if (e2 < dx) { err += dx; py += sy; }
    }

    int rowSize = ((bmp.infoHeader.bitCount * bmp.infoHeader.width + 31) / 32) * 4;
    auto worker = [&](int start, int end) {
        for (int i = start; i < end; ++i) {
            int x = points[i].first;
            int y = points[i].second;
            for (int tx = -thickness / 2; tx <= thickness / 2; ++tx) {
                for (int ty = -thickness / 2; ty <= thickness / 2; ++ty) {
                    int px = x + tx;
                    int py = y + ty;
                    if (px >= 0 && px < bmp.infoHeader.width && py >= 0 && py < bmp.infoHeader.height) {
                        int index = py * rowSize + px * 3;
                        bmp.data[index] = b;
                        bmp.data[index + 1] = g;
                        bmp.data[index + 2] = r;
                    }
                }
            }
        }
        };

    std::vector<std::thread> threads;
    int chunk = points.size() / numThreads;
    for (int i = 0; i < numThreads; ++i) {
        int start = i * chunk;
        int end = (i == numThreads - 1) ? points.size() : start + chunk;
        threads.emplace_back(worker, start, end);
    }

    for (auto& t : threads) t.join();
}

std::vector<Point> readPointsFromCSV(const std::string& filename) {
    std::vector<Point> points;
    std::ifstream file(filename);

    if (!file.is_open()) {
        throw std::runtime_error("Unable to open CSV file.");
    }

    std::string line;
    while (std::getline(file, line)) {
        std::stringstream ss(line);
        int x, y;
        char comma;
        ss >> x >> comma >> y;

        points.push_back({ x, y });
    }

    return points;
}

std::vector<Point> generateRandomPoints(int numPoints, int maxX, int maxY) {
    std::vector<Point> points;
    std::srand(static_cast<unsigned>(std::time(0)));

    for (int i = 0; i < numPoints; ++i) {
        points.push_back({ std::rand() % maxX, std::rand() % maxY });
    }

    return points;
}
void drawPolyline(BMP& bmp, const std::vector<Point>& points, int thickness, uint8_t r, uint8_t g, uint8_t b) {
    for (size_t i = 0; i < points.size() - 1; ++i) {
        drawLine(bmp, points[i].x, points[i].y, points[i + 1].x, points[i + 1].y, thickness, r, g, b);
    }
}

void drawPolylineThreaded(BMP& bmp, const std::vector<Point>& points, int thickness,
    uint8_t r, uint8_t g, uint8_t b, int numThreads) {
    int numSegments = points.size() - 1;
    if (numSegments <= 0 || numThreads <= 0) return;

    std::thread* threads = new std::thread[numThreads];  // dynamiczna tablica

    int segmentsPerThread = numSegments / numThreads;
    int remaining = numSegments % numThreads;

    int start = 0;
    for (int t = 0; t < numThreads; ++t) {
        int end = start + segmentsPerThread + (t < remaining ? 1 : 0);

        threads[t] = std::thread([=, &bmp, &points]() {
            for (int i = start; i < end; ++i) {
                drawLine(bmp, points[i].x, points[i].y,
                    points[i + 1].x, points[i + 1].y,
                    thickness, r, g, b);
            }
            });

        start = end;
    }

    for (int t = 0; t < numThreads; ++t) {
        if (threads[t].joinable()) threads[t].join();
    }

    delete[] threads;  // zwolnij pamięć
}

void drawPolylineOpenMP(BMP& bmp, const std::vector<Point>& points, int thickness,
    uint8_t r, uint8_t g, uint8_t b, int numThreads) {
    int numSegments = points.size() - 1;

   // omp_set_num_threads(numThreads);

#pragma omp parallel for
    for (int i = 0; i < numSegments; ++i) {
        drawLine(bmp, points[i].x, points[i].y,
            points[i + 1].x, points[i + 1].y,
            thickness, r, g, b);
    }
}
// obliczenia automatyczne
void test_line_rasterization_lengths() {
    const int imgSize = 512;
    const char outCSV[] = "line.csv";
    const std::vector<int> lengths = { 300 };
    const int minThickness = 0;
    const int maxThickness = 300;
    const int step = 25;
    std::ofstream results(outCSV);
    results << "thickness,sequential\n";

    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<> dist(0, imgSize - 1);

    for (int len : lengths) {
        for (int thickness = minThickness; thickness <= maxThickness; thickness+=step) {
            // Szukaj pary punktów o zadanej długości (z tolerancją)
            int x1, y1, x2, y2;
            double foundLength = 0.0;
            int attempts = 0;
            do {
                x1 = dist(gen);
                y1 = dist(gen);
                x2 = dist(gen);
                y2 = dist(gen);
                foundLength = std::sqrt((x2 - x1) * (x2 - x1) + (y2 - y1) * (y2 - y1));
                ++attempts;
                // Zabezpieczenie przed nieskończoną pętlą
                if (attempts > 10000) break;
            } while (std::abs(foundLength - len) > 0.5);

            // Jeśli nie znaleziono odpowiedniej pary, pomiń ten przypadek
            if (std::abs(foundLength - len) > 0.5)
                continue;

            BMP bmp;
            bmp.infoHeader.width = imgSize;
            bmp.infoHeader.height = imgSize;
            bmp.infoHeader.bitCount = 24;
            int rowSize = ((bmp.infoHeader.bitCount * bmp.infoHeader.width + 31) / 32) * 4;
            int dataSize = rowSize * imgSize;
            bmp.data = new uint8_t[dataSize];
            std::memset(bmp.data, 255, dataSize); // white background
            int t;
			if (thickness <= 0) {
                t = 1;
			}
			t = thickness;
            LARGE_INTEGER start, end, frequency;
            QueryPerformanceFrequency(&frequency);
            QueryPerformanceCounter(&start);

            drawLine(bmp, x1, y1, x2, y2, t, 255, 0, 0);

            QueryPerformanceCounter(&end);
            double t_seq = (end.QuadPart - start.QuadPart) * 1000.0 / frequency.QuadPart;

            results << t << "," << t_seq << "\n";

            delete[] bmp.data;
        }
    }

    results.close();
    std::cout << "Wyniki zapisane do " << outCSV << "\n";
}
void test_line_rasterization_lengths_threaded() {
    const int imgSize = 512;
    const char outCSV[] = "line_threads.csv";
    const int length = 300;
    const int minThickness = 1;
    const int maxThickness = 300;
    const int step = 25;
    const std::vector<int> threadCounts = { 2, 4, 7, 8, 9, 10, 16, 32, 64, 128, 256 };

    std::ofstream results(outCSV);
    results << "thickness,threads,threaded\n";

    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<> dist(0, imgSize - 1);

    for (int thickness = 1, i = 0; thickness <= maxThickness; thickness = (i == 0 ? 25 : thickness + step), ++i) {
        for (int numThreads : threadCounts) {
            // Find a pair of points with the required length (tolerance ±0.5)
            int x1, y1, x2, y2;
            double foundLength = 0.0;
            int attempts = 0;
            do {
                x1 = dist(gen);
                y1 = dist(gen);
                x2 = dist(gen);
                y2 = dist(gen);
                foundLength = std::sqrt((x2 - x1) * (x2 - x1) + (y2 - y1) * (y2 - y1));
                ++attempts;
                if (attempts > 10000) break;
            } while (std::abs(foundLength - length) > 0.5);

            if (std::abs(foundLength - length) > 0.5)
                continue;

            BMP bmp;
            bmp.infoHeader.width = imgSize;
            bmp.infoHeader.height = imgSize;
            bmp.infoHeader.bitCount = 24;
            int rowSize = ((bmp.infoHeader.bitCount * bmp.infoHeader.width + 31) / 32) * 4;
            int dataSize = rowSize * imgSize;
            bmp.data = new uint8_t[dataSize];
            std::memset(bmp.data, 255, dataSize); // white background

            LARGE_INTEGER start, end, frequency;
            QueryPerformanceFrequency(&frequency);
            QueryPerformanceCounter(&start);

            drawLineThreaded(bmp, x1, y1, x2, y2, thickness, 255, 0, 0, numThreads);

            QueryPerformanceCounter(&end);
            double t_thr = (end.QuadPart - start.QuadPart) * 1000.0 / frequency.QuadPart;

            results << thickness << "," << numThreads << "," << t_thr << "\n";

            delete[] bmp.data;
        }
    }

    results.close();
    std::cout << "Wyniki (threaded) zapisane do " << outCSV << "\n";
}

void test_line_rasterization_lengths_openmp() {
    const int imgSize = 512;
    const char outCSV[] = "line_openmp.csv";
    const int length = 300;
    const int minThickness = 1;
    const int maxThickness = 300;
    const int step = 25;
    const std::vector<int> threadCounts = { 2, 4, 7, 8, 9, 10, 16, 32, 64, 128, 256 };

    std::ofstream results(outCSV);
    results << "thickness,threads,openmp\n";

    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<> dist(0, imgSize - 1);

    for (int thickness = 1, i = 0; thickness <= maxThickness; thickness = (i == 0 ? 25 : thickness + step), ++i) {
        for (int numThreads : threadCounts) {
            // Find a pair of points with the required length (tolerance ±0.5)
            int x1, y1, x2, y2;
            double foundLength = 0.0;
            int attempts = 0;
            do {
                x1 = dist(gen);
                y1 = dist(gen);
                x2 = dist(gen);
                y2 = dist(gen);
                foundLength = std::sqrt((x2 - x1) * (x2 - x1) + (y2 - y1) * (y2 - y1));
                ++attempts;
                if (attempts > 10000) break;
            } while (std::abs(foundLength - length) > 0.5);

            if (std::abs(foundLength - length) > 0.5)
                continue;

            BMP bmp;
            bmp.infoHeader.width = imgSize;
            bmp.infoHeader.height = imgSize;
            bmp.infoHeader.bitCount = 24;
            int rowSize = ((bmp.infoHeader.bitCount * bmp.infoHeader.width + 31) / 32) * 4;
            int dataSize = rowSize * imgSize;
            bmp.data = new uint8_t[dataSize];
            std::memset(bmp.data, 255, dataSize); // white background

            LARGE_INTEGER start, end, frequency;
            QueryPerformanceFrequency(&frequency);
            QueryPerformanceCounter(&start);

            drawLineOpenMP(bmp, x1, y1, x2, y2, thickness, 255, 0, 0, numThreads);

            QueryPerformanceCounter(&end);
            double t_omp = (end.QuadPart - start.QuadPart) * 1000.0 / frequency.QuadPart;

            results << thickness << "," << numThreads << "," << t_omp << "\n";

            delete[] bmp.data;
        }
    }

    results.close();
    std::cout << "Wyniki (openmp) zapisane do " << outCSV << "\n";
}

void PolyLineSeq() {
    try {
        // Read points from CSV
        std::vector<Point> points = readPointsFromCSV("punkty.csv");

        // Prepare CSV output
        std::ofstream results("poly.csv");

        // Loop over thickness values
        for (double thick = 1.0; thick <= 50; thick += 1) {
            BMP bmp1 = readBMP("input.bmp");

            LARGE_INTEGER start, end, frequency;
            QueryPerformanceFrequency(&frequency);
            QueryPerformanceCounter(&start);

            drawPolyline(bmp1, points, static_cast<int>(thick), 255, 0, 0);

            QueryPerformanceCounter(&end);
            double t_seq = (end.QuadPart - start.QuadPart) * 1000.0 / frequency.QuadPart;

            // Save output image (optional, can be commented out if not needed)
            // saveBMP(outSeq, bmp1);
            delete[] bmp1.data;

            // Write to CSV: thickness, sequential, threaded, openmp
            results << thick << "," << t_seq << "\n";
        }

        results.close();
        std::cout << "Sequential timing results saved to " << "poly.csv" << "\n";
    }
    catch (const std::exception& e) {
        std::cerr << "Blad: " << e.what() << "\n";
    }
}
void PolyLineThreaded(const std::vector<Point>& points) {
    const char* outCSV = "poly_threaded.csv";
    const std::vector<int> threadCounts = { 2, 4, 7, 8, 9, 10, 16, 32, 64, 128, 256 };

    std::ofstream results(outCSV);
    results << "thickness,threads,threaded\n";

    for (int thick = 1; thick <= 50; ++thick) {
        for (int numThreads : threadCounts) {
            BMP bmp = readBMP("input.bmp");

            LARGE_INTEGER start, end, frequency;
            QueryPerformanceFrequency(&frequency);
            QueryPerformanceCounter(&start);

            drawPolylineThreaded(bmp, points, thick, 0, 255, 0, numThreads);

            QueryPerformanceCounter(&end);
            double t_thr = (end.QuadPart - start.QuadPart) * 1000.0 / frequency.QuadPart;

            delete[] bmp.data;

            results << thick << "," << numThreads << "," << t_thr << "\n";
        }
    }

    results.close();
    std::cout << "Threaded timing results saved to " << outCSV << "\n";
}

void PolyLineOpenMP(const std::vector<Point>& points) {
    const char* outCSV = "poly_openmp.csv";
    const std::vector<int> threadCounts = { 2, 4, 7, 8, 9, 10, 16, 32, 64, 128, 256 };

    std::ofstream results(outCSV);
    results << "thickness,threads,openmp\n";

    for (int thick = 1; thick <= 50; ++thick) {
        for (int numThreads : threadCounts) {
            BMP bmp = readBMP("input.bmp");

            LARGE_INTEGER start, end, frequency;
            QueryPerformanceFrequency(&frequency);
            QueryPerformanceCounter(&start);

            drawPolylineOpenMP(bmp, points, thick, 0, 0, 255, numThreads);

            QueryPerformanceCounter(&end);
            double t_omp = (end.QuadPart - start.QuadPart) * 1000.0 / frequency.QuadPart;

            delete[] bmp.data;

            results << thick << "," << numThreads << "," << t_omp << "\n";
        }
    }

    results.close();
    std::cout << "OpenMP timing results saved to " << outCSV << "\n";
}

int main() {
    const char filename[] = "test2.bmp";
    const char outSeq[] = "out_seq.bmp";
    const char outThreaded[] = "out_threaded.bmp";
    const char outOmp[] = "out_omp.bmp";
    int x1, y1, x2, y2, thick;
    int numPoints, maxX, maxY;
    std::string csvFile;
    int numThreads;
    // automating changes
    std::vector<Point> points;
    points = readPointsFromCSV("punkty.csv");
    const char resultsFile[] = "results.csv";


    //test_line_rasterization_lengths();
    //test_line_rasterization_lengths_threaded();
    //test_line_rasterization_lengths_openmp();

    //PolyLineSeq(); // sequential
    PolyLineThreaded(points); // threaded
    PolyLineOpenMP(points);
    

/*
    try {
        BMP bmp = readBMP(filename);
        maxX = bmp.infoHeader.width;
        maxY = bmp.infoHeader.height;

        std::cout << "Wybierz scenariusz:" << std::endl;
        std::cout << "1. Wczytanie punktów z pliku CSV" << std::endl;
        std::cout << "2. Generowanie punktow losowo" << std::endl;
        int choice;
        std::cin >> choice;

        std::vector<Point> points;

        if (choice == 1) {
            std::cout << "Podaj nazwe pliku CSV: ";
            std::cin >> csvFile;
            points = readPointsFromCSV(csvFile);
        }
        else if (choice == 2) {
            std::cout << "Podaj liczbe punktow: ";
            std::cin >> numPoints;
            points = generateRandomPoints(numPoints, maxX, maxY);
        }
        else {
            throw std::invalid_argument("Niepoprawny wybor.");
        }

        std::cout << "Podaj grubosc lamanej [px]: ";
        std::cin >> thick;

        std::cout << "Podaj liczbe grup: ";
        std::cin >> numThreads;

        // Rysowanie łamanej w różnych wersjach

        // SEKWENCYJNIE
        BMP bmp1 = readBMP(filename);
        LARGE_INTEGER start, end, frequency;
        QueryPerformanceFrequency(&frequency);
        QueryPerformanceCounter(&start);
        drawPolyline(bmp1, points, thick, 255, 0, 0);
        QueryPerformanceCounter(&end);
        double t_seq = (end.QuadPart - start.QuadPart) * 1000.0 / frequency.QuadPart;
        std::cout << "Czas SEKWENCYJNY: " << t_seq << " ms\n";
        saveBMP(outSeq, bmp1);
        delete[] bmp1.data;

        // THREADED
        BMP bmp2 = readBMP(filename);
        QueryPerformanceCounter(&start);
        drawPolylineThreaded(bmp2, points, thick, 0, 255, 0, numThreads);
        QueryPerformanceCounter(&end);
        double t_thr = (end.QuadPart - start.QuadPart) * 1000.0 / frequency.QuadPart;
        std::cout << "Czas THREADED(" << numThreads << " watki): " << t_thr << " ms\n";
        saveBMP(outThreaded, bmp2);
        delete[] bmp2.data;

        // OPENMP
        BMP bmp3 = readBMP(filename);
        QueryPerformanceCounter(&start);
        drawPolylineOpenMP(bmp3, points, thick, 0, 0, 255, numThreads);
        QueryPerformanceCounter(&end);
        double t_omp = (end.QuadPart - start.QuadPart) * 1000.0 / frequency.QuadPart;
        std::cout << "Czas OpenMP: " << t_omp << " ms\n";
        saveBMP(outOmp, bmp3);
        delete[] bmp3.data;

    }
    catch (const std::exception& e) {
        std::cerr << "Blad: " << e.what() << "\n";
    }
    */
    return 0;
}

// problemem jest granulacja i obliczenie czy powtarzalnosc testu wplywa na czasy.