#include <iostream>
#include <cmath>
#include <limits>
#include <vector>
#include <fstream>
#include <algorithm>

struct Vector3 {
    float x, y, z;

    Vector3 operator+(const Vector3 &ref) const {
        return {x + ref.x, y + ref.y, z + ref.z};
    }
    Vector3 operator-(const Vector3 &ref) const {
        return {x - ref.x, y - ref.y, z - ref.z};
    }
    Vector3 operator*(float s) const {
        return {x * s, y * s, z * s};
    }

    Vector3 operator*(const Vector3 &ref) const {
        return {x * ref.x, y * ref.y, z * ref.z};
    }

    Vector3 operator/(float s) const {
        if (s == 0)
            return {0, 0, 0};
        return {x / s, y / s, z / s};
    }

    float dot(const Vector3 &ref) const {
        return x * ref.x + y * ref.y + z * ref.z;
    }

    Vector3 cross(const Vector3 &ref) const {
        return {
            y * ref.z - z * ref.y,
            z * ref.x - x * ref.z,
            x * ref.y - y * ref.x};
    }
    float length() const {
        return std::sqrt(x * x + y * y + z * z);
    }
    Vector3 normalized() const {
        float len = length();
        if (len < std::numeric_limits<float>::epsilon())
            return {0, 0, 0};
        return *this / len;
    }
};

std::ostream &operator<<(std::ostream &os, const Vector3 &p) {
    os << "(" << p.x << ", " << p.y << ", " << p.z << ")";
    return os;
}

struct light {
    Vector3 position;
    Vector3 color;
    float intensity;
};

struct Ray {
    Vector3 origin;
    Vector3 direction;
    Ray(const Vector3 &o, const Vector3 &d)
        : origin(o), direction(d.normalized()) {}
};

struct Sphere {
    Vector3 center;
    float radius;
    Vector3 color;
    float k_s;
    float shinignes;

    bool intersect(const Ray &ray, float &t0) const {
        Vector3 L = center - ray.origin;
        float tca = L.dot(ray.direction);
        if (tca < 0)
            return false;
        float d2 = L.dot(L) - tca * tca;
        float r2 = radius * radius;
        if (d2 > r2)
            return false;
        float thc = std::sqrt(r2 - d2);
        float t1 = tca + thc;
        t0 = tca - thc;
        if (t0 < 0)
            t0 = t1;
        if (t0 < 0)
            return false;
        return true;
    }
};

Vector3 ray_cast(std::vector<Sphere> &spheres, float &nearest_t, Ray &ray, Vector3 &pixel) {
    std::vector<light> lights = {
        {{80, 100, 50}, {1.0, 0.95, 0.8}, 2.5},
        {{-60, 40, 20}, {0.6, 0.8, 1.0}, 1.0},
        {{0, 50, 100}, {1, 1, 1}, 0.5},
        {{0, -100, -10}, {1.0, 0.8, 0.6}, 0.2}};

    Vector3 total_color = {0, 0, 0};

    for (const auto &sphere : spheres) {
        float t0;
        if (sphere.intersect(ray, t0)) {
            if (t0 < nearest_t) {
                for (auto &swiatlo : lights) {
                    nearest_t = t0;
                    Vector3 P = ray.origin + (ray.direction * t0);
                    Vector3 N = (P - sphere.center).normalized();
                    Vector3 L = swiatlo.position - P;
                    Vector3 V = (ray.origin - P).normalized();
                    L = L.normalized();
                    float distance = L.length();
                    N = N.normalized();
                    Vector3 R = N * 2 * N.dot(L) - L;
                    float I = std::max(0.0f, N.dot(L));
                    I = I / (distance * distance);
                    I *= swiatlo.intensity;
                    float specular = std::pow((std::max(0.0f, R.dot(V))), sphere.shinignes) * sphere.k_s;

                    total_color = total_color + ((sphere.color * swiatlo.color * I) + (swiatlo.color * specular));
                }

                return total_color;
            }
        }
    }
    return pixel;
}

void render() {
    const int width = 1920;
    const int height = 1080;
    Vector3 camera{0, 0, 0};
    float fov = M_PI / 2.0f;
    std::vector<Vector3> framebuffer(width * height);
    std::vector<Sphere> spheres = {
        {{2, 0, -3}, 1, {0.67, 0, 0}, 0.12, 30},
        {{-5, -4, -6}, 1.4, {0.67, 0, 0}, 0.35, 120},
        {{-1, 3, -30}, 10, {0.1, 0.41, 0.1}, 0.95, 350},
        {{-3, 3, -6}, 3, {0.8, 0.2, 0.2}, 1, 800}

    };

    for (int j = 0; j < height; ++j) {
        for (int i = 0; i < width; ++i) {
            float x = (2 * (i + 0.5f) / float(width) - 1) * std::tan(fov / 2.0f) * width / float(height);
            float y = -(2 * (j + 0.5f) / float(height) - 1) * std::tan(fov / 2.0f);
            Ray ray(camera, Vector3{x, y, -1});
            Vector3 pixel_color = {0.1f, 0.1f, 0.1f};
            float nearest_t = std::numeric_limits<float>::max();

            pixel_color = ray_cast(spheres, nearest_t, ray, pixel_color);

            framebuffer[i + j * width] = pixel_color;
        }
    }

    std::ofstream ofs("out.ppm", std::ios::binary);
    ofs << "P6\n"
        << width << " " << height << "\n255\n";
    for (auto &c : framebuffer) {
        ofs << (unsigned char)(255 * std::clamp(c.x, 0.0f, 1.0f));
        ofs << (unsigned char)(255 * std::clamp(c.y, 0.0f, 1.0f));
        ofs << (unsigned char)(255 * std::clamp(c.z, 0.0f, 1.0f));
    }
    ofs.close();
}

int main() {
    render();
    return 0;
}
