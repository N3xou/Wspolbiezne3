#include <iostream>
#include <fstream>
#include <vector>
#include <queue>
#include <cstdlib>
#include <ctime>
#include <windows.h>
#pragma pack(push, 1)
#include <thread>
#include <omp.h>
#include <mutex>
#include <atomic>
struct BMPHeader {
    uint16_t type;
    uint32_t size;
    uint16_t reserved1, reserved2;
    uint32_t offset;
};

struct DIBHeader {
    uint32_t size;
    int32_t width, height;
    uint16_t planes;
    uint16_t bitCount;
    uint32_t compression;
    uint32_t imageSize;
    int32_t xPelsPerMeter, yPelsPerMeter;
    uint32_t clrUsed, clrImportant;
};

#pragma pack(pop)

struct RGB {
    uint8_t b, g, r;
    bool operator==(const RGB& other) const {
        return r == other.r && g == other.g && b == other.b;
    }
    bool operator!=(const RGB& other) const {
        return !(*this == other);
    }
};

const RGB BLACK = { 0, 0, 0 };
const RGB WHITE = { 255, 255, 255 };

std::vector<std::vector<RGB>> image;
std::vector<std::vector<bool>> visited;
int width, height;

std::mutex visited_mutex;

std::vector<std::pair<int, int>> get_starting_points() {
    std::vector<std::pair<int, int>> points;
    for (int y = 0; y < height; ++y)
        for (int x = 0; x < width; ++x)
            if (!visited[y][x] && image[y][x] != BLACK)
                points.emplace_back(x, y);
    return points;
}
bool is_inside(int x, int y) {
    return x >= 0 && y >= 0 && x < width && y < height;
}

bool is_closed_shape(int x, int y) {
    std::queue<std::pair<int, int>> q;
    std::vector<std::pair<int, int>> points;
    bool touches_edge = false;

    q.push({ x, y });
    visited[y][x] = true;
    RGB fill_color = BLACK;
    for (int dx = -1; dx <= 1; ++dx) {
        for (int dy = -1; dy <= 1; ++dy) {
            int nx = x + dx;
            int ny = y + dy;
            if ((dx != 0 || dy != 0) && is_inside(nx, ny)) {
                if (image[ny][nx] != WHITE) { // Assuming WHITE is the background
                    fill_color = image[ny][nx];
                    break;
                }
            }
        }
    }
    while (!q.empty()) {
        std::pair<int, int> front = q.front(); q.pop();
        int cx = front.first;
        int cy = front.second;
        points.push_back({ cx, cy });

        if (cx == 0 || cy == 0 || cx == width - 1 || cy == height - 1) {
            touches_edge = true;
        }

        for (int dx : {-1, 1, 0, 0}) {
            for (int dy : {0, 0, -1, 1}) {
                int nx = cx + dx;
                int ny = cy + dy;

                if (is_inside(nx, ny) && !visited[ny][nx]) {
                    RGB color = image[ny][nx];
                    if (color.r == 255 && color.g == 255 && color.b == 255) {
                        visited[ny][nx] = true;
                        q.push({ nx, ny });
                    }
                }
            }
        }
    }

    if (!touches_edge) {
        // Flood-fill with a random color


        for (const auto& point : points) {
            int px = point.first;
            int py = point.second;
            image[py][px] = fill_color;
        }
    }

    return !touches_edge;
}
void fill_shapes_threaded(int num_threads) {
    auto points = get_starting_points();
    size_t idx = 0;

    auto worker = [&]() {
        while (true) {
            size_t i;
            {
                //std::lock_guard<std::mutex> lock(visited_mutex);
                if (idx >= points.size()) return;
                i = idx++;
            }
            int x = points[i].first;
            int y = points[i].second;
            bool should_fill = false;
            {
                //std::lock_guard<std::mutex> lock(visited_mutex);
                if (!visited[y][x] && image[y][x] != BLACK) {
                    visited[y][x] = true;
                    should_fill = true;
                }
            }
            if (should_fill) {
                is_closed_shape(x, y);
            }
        }
        };

    std::vector<std::thread> threads;
    for (int t = 0; t < num_threads; ++t)
        threads.emplace_back(worker);
    for (auto& th : threads)
        th.join();
}



void fill_shapes_openmp(int num_threads) {
    auto points = get_starting_points();
    size_t n = points.size();

#pragma omp parallel for num_threads(num_threads) schedule(dynamic)
    for (int i = 0; i < static_cast<int>(n); ++i) {
        int x = points[i].first;
        int y = points[i].second;
        bool should_fill = false;
        {
            //std::lock_guard<std::mutex> lock(visited_mutex);
            if (!visited[y][x] && image[y][x] != BLACK) {
                visited[y][x] = true;
                should_fill = true;
            }
        }
        if (should_fill) {
            is_closed_shape(x, y);
        }
    }
}





void read_bmp(const std::string& filename) {
    std::ifstream in(filename, std::ios::binary);
    if (!in) {
        std::cerr << "Cannot open file\n";
        exit(1);
    }

    BMPHeader bmp;
    DIBHeader dib;

    in.read(reinterpret_cast<char*>(&bmp), sizeof(bmp));
    in.read(reinterpret_cast<char*>(&dib), sizeof(dib));

    if (bmp.type != 0x4D42 || dib.bitCount != 24) {
        std::cerr << "Unsupported BMP format\n";
        exit(1);
    }

    width = dib.width;
    height = dib.height;

    image.resize(height, std::vector<RGB>(width));
    visited.resize(height, std::vector<bool>(width, false));

    int row_padded = (width * 3 + 3) & (~3);
    std::vector<uint8_t> row(row_padded);

    in.seekg(bmp.offset, std::ios::beg);

    for (int y = height - 1; y >= 0; --y) {
        in.read(reinterpret_cast<char*>(row.data()), row_padded);
        for (int x = 0; x < width; ++x) {
            image[y][x] = { row[x * 3], row[x * 3 + 1], row[x * 3 + 2] };
        }
    }

    in.close();
}

void write_bmp(const std::string& filename) {
    std::ofstream out(filename, std::ios::binary);
    if (!out) {
        std::cerr << "Cannot write file\n";
        exit(1);
    }

    BMPHeader bmp = { 0x4D42, 0, 0, 0, 54 };
    DIBHeader dib = { 40, width, height, 1, 24, 0, 0, 2835, 2835, 0, 0 };

    int row_padded = (width * 3 + 3) & (~3);
    bmp.size = bmp.offset + row_padded * height;

    out.write(reinterpret_cast<char*>(&bmp), sizeof(bmp));
    out.write(reinterpret_cast<char*>(&dib), sizeof(dib));

    std::vector<uint8_t> row(row_padded, 0);

    for (int y = height - 1; y >= 0; --y) {
        for (int x = 0; x < width; ++x) {
            RGB pix = image[y][x];
            row[x * 3] = pix.b;
            row[x * 3 + 1] = pix.g;
            row[x * 3 + 2] = pix.r;
        }
        out.write(reinterpret_cast<char*>(row.data()), row_padded);
    }

    out.close();
}

int main() {
    LARGE_INTEGER start, end, frequency;
    int num_parts;
    QueryPerformanceFrequency(&frequency);
    srand(time(0));
    read_bmp("ksztalty.bmp");
 
    std::cout << "Podaj liczbe czesci/watkow (1-500): ";
    std::cin >> num_parts;
    QueryPerformanceCounter(&start);
    for (int y = 0; y < height; ++y)
        for (int x = 0; x < width; ++x)
            if (!visited[y][x] && image[y][x] != BLACK)
                is_closed_shape(x, y);
    QueryPerformanceCounter(&end);
    double t_seq = (end.QuadPart - start.QuadPart) * 1000.0 / frequency.QuadPart;
    std::cout << "Czas SEKWENCYJNY: " << t_seq << " ms\n";
    std::cout << "Czas SEKWENCYJNY: " << (t_seq / 1000.0) << " s\n";
    write_bmp("output.bmp");
    std::cout << "Saved to " << "output.bmp" << "\n";

    // --- THREADED ---
    // Reset visited and reload image
    read_bmp("ksztalty.bmp");
    for (auto& row : visited)
        std::fill(row.begin(), row.end(), false);

    QueryPerformanceCounter(&start);
    fill_shapes_threaded(num_parts);
    QueryPerformanceCounter(&end);
    double t_threaded = (end.QuadPart - start.QuadPart) * 1000.0 / frequency.QuadPart;
    std::cout << "Czas WATKI: " << t_threaded << " ms\n";
    std::cout << "Czas WATKI: " << (t_threaded / 1000.0) << " s\n";
    write_bmp("outputthreaded.bmp");
    std::cout << "Saved to outputthreaded.bmp\n";

    // --- OPENMP ---
    // Reset visited and reload image
    read_bmp("ksztalty.bmp");
    for (auto& row : visited)
        std::fill(row.begin(), row.end(), false);

    QueryPerformanceCounter(&start);
    fill_shapes_openmp(num_parts);
    QueryPerformanceCounter(&end);
    double t_openmp = (end.QuadPart - start.QuadPart) * 1000.0 / frequency.QuadPart;
    std::cout << "Czas OPENMP: " << t_openmp << " ms\n";
    std::cout << "Czas OPENMP: " << (t_openmp / 1000.0) << " s\n";
    write_bmp("outputopenmp.bmp");
    std::cout << "Saved to outputopenmp.bmp\n";


    return 0;
}
