#include <iostream>
#include <fstream>
#include <vector>
#include <queue>
#include <cstdlib>
#include <ctime>

#pragma pack(push, 1)

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

bool is_inside(int x, int y) {
    return x >= 0 && y >= 0 && x < width && y < height;
}

bool is_closed_shape(int x, int y, const RGB& border_color) {
    std::queue<std::pair<int, int>> q;
    std::vector<std::pair<int, int>> points;
    bool touches_edge = false;

    q.push({ x, y });
    visited[y][x] = true;

    while (!q.empty()) {
        auto [cx, cy] = q.front(); q.pop();
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
                    if (color != border_color) {
                        visited[ny][nx] = true;
                        q.push({ nx, ny });
                    }
                }
            }
        }
    }

    if (!touches_edge) {
        // Flood-fill with a random color
        RGB fill_color = {
            static_cast<uint8_t>(rand() % 200 + 30),
            static_cast<uint8_t>(rand() % 200 + 30),
            static_cast<uint8_t>(rand() % 200 + 30)
        };

        for (auto [px, py] : points) {
            image[py][px] = fill_color;
        }
    }

    return !touches_edge;
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

int main(int argc, char* argv[]) {
    srand(time(0));

    if (argc != 3) {
        std::cout << "Usage: ./fill_shapes input.bmp output.bmp\n";
        return 1;
    }

    read_bmp(argv[1]);

    for (int y = 0; y < height; ++y)
        for (int x = 0; x < width; ++x)
            if (!visited[y][x] && image[y][x] != BLACK)
                is_closed_shape(x, y, BLACK);

    write_bmp(argv[2]);
    std::cout << "Saved to " << argv[2] << "\n";

    return 0;
}
