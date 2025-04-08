#include <iostream>
#include <cmath>
#include <utility>
#include <vector>
#include <fstream>
#include <tuple>

struct CelestialBody
{
    std::string name;
    double mass;
    std::tuple<double, double, double> position;  // x, y, z
    std::tuple<double, double, double> velocity;  // vx, vy, vz
    std::tuple<double, double, double> acceleration;  // ax, ay, az
};


struct TrajectoryPoint
{
    long double time;
    std::string name;
    std::tuple<double, double, double> pos;
};


double distance(std::tuple<double, double, double> pos1, std::tuple<double, double, double> pos2)
{
    double dx = std::get<0>(pos2) - std::get<0>(pos1);
    double dy = std::get<1>(pos2) - std::get<1>(pos1);
    double dz = std::get<2>(pos2) - std::get<2>(pos1);
    return sqrt(dx * dx + dy * dy + dz * dz);
}


void calculate_accelerations(std::vector<CelestialBody>& bodies)
{
    const double g = 6.6743e-11;

    // ускорения
    for (auto& body : bodies)
    {
        body.acceleration = std::make_tuple(0.0, 0.0, 0.0);
    }

    // взаимодействия между всеми парами тел
    for (size_t i = 0; i < bodies.size(); ++i)
    {
        for (size_t j = i + 1; j < bodies.size(); ++j)
        {
            double r = distance(bodies[i].position, bodies[j].position);
            double F = g * bodies[i].mass * bodies[j].mass / (r * r * r);

            double dx = std::get<0>(bodies[j].position) - std::get<0>(bodies[i].position);
            double dy = std::get<1>(bodies[j].position) - std::get<1>(bodies[i].position);
            double dz = std::get<2>(bodies[j].position) - std::get<2>(bodies[i].position);

            std::get<0>(bodies[i].acceleration) += F * dx;
            std::get<1>(bodies[i].acceleration) += F * dy;
            std::get<2>(bodies[i].acceleration) += F * dz;

            std::get<0>(bodies[j].acceleration) -= F * dx;
            std::get<1>(bodies[j].acceleration) -= F * dy;
            std::get<2>(bodies[j].acceleration) -= F * dz;
        }
    }

    for (auto& body : bodies)
    {
        std::get<0>(body.acceleration) /= body.mass;
        std::get<1>(body.acceleration) /= body.mass;
        std::get<2>(body.acceleration) /= body.mass;
    }
}


void runge_kutta_step(std::vector<CelestialBody>& bodies, double dt)
{
    struct State
    {
        std::tuple<double, double, double> pos;
        std::tuple<double, double, double> vel;
    };

    std::vector<State> k1(bodies.size()), k2(bodies.size()), k3(bodies.size()), k4(bodies.size());
    std::vector<CelestialBody> temp_bodies = bodies;

    calculate_accelerations(temp_bodies);

    for (size_t i = 0; i < bodies.size(); ++i)
    {
        k1[i].pos = std::make_tuple(
            std::get<0>(temp_bodies[i].velocity) * dt,
            std::get<1>(temp_bodies[i].velocity) * dt,
            std::get<2>(temp_bodies[i].velocity) * dt
        );
        k1[i].vel = std::make_tuple(
            std::get<0>(temp_bodies[i].acceleration) * dt,
            std::get<1>(temp_bodies[i].acceleration) * dt,
            std::get<2>(temp_bodies[i].acceleration) * dt
        );

        std::get<0>(temp_bodies[i].position) += std::get<0>(k1[i].pos) / 2;
        std::get<1>(temp_bodies[i].position) += std::get<1>(k1[i].pos) / 2;
        std::get<2>(temp_bodies[i].position) += std::get<2>(k1[i].pos) / 2;
        std::get<0>(temp_bodies[i].velocity) += std::get<0>(k1[i].vel) / 2;
        std::get<1>(temp_bodies[i].velocity) += std::get<1>(k1[i].vel) / 2;
        std::get<2>(temp_bodies[i].velocity) += std::get<2>(k1[i].vel) / 2;
    }

    calculate_accelerations(temp_bodies);

    for (size_t i = 0; i < bodies.size(); ++i)
    {
        k2[i].pos = std::make_tuple(
            std::get<0>(temp_bodies[i].velocity) * dt,
            std::get<1>(temp_bodies[i].velocity) * dt,
            std::get<2>(temp_bodies[i].velocity) * dt
        );
        k2[i].vel = std::make_tuple(
            std::get<0>(temp_bodies[i].acceleration) * dt,
            std::get<1>(temp_bodies[i].acceleration) * dt,
            std::get<2>(temp_bodies[i].acceleration) * dt
        );

        std::get<0>(temp_bodies[i].position) = std::get<0>(bodies[i].position) + std::get<0>(k2[i].pos) / 2;
        std::get<1>(temp_bodies[i].position) = std::get<1>(bodies[i].position) + std::get<1>(k2[i].pos) / 2;
        std::get<2>(temp_bodies[i].position) = std::get<2>(bodies[i].position) + std::get<2>(k2[i].pos) / 2;
        std::get<0>(temp_bodies[i].velocity) = std::get<0>(bodies[i].velocity) + std::get<0>(k2[i].vel) / 2;
        std::get<1>(temp_bodies[i].velocity) = std::get<1>(bodies[i].velocity) + std::get<1>(k2[i].vel) / 2;
        std::get<2>(temp_bodies[i].velocity) = std::get<2>(bodies[i].velocity) + std::get<2>(k2[i].vel) / 2;
    }


    calculate_accelerations(temp_bodies);

    for (size_t i = 0; i < bodies.size(); ++i)
    {
        k3[i].pos = std::make_tuple(
            std::get<0>(temp_bodies[i].velocity) * dt,
            std::get<1>(temp_bodies[i].velocity) * dt,
            std::get<2>(temp_bodies[i].velocity) * dt
        );
        k3[i].vel = std::make_tuple(
            std::get<0>(temp_bodies[i].acceleration) * dt,
            std::get<1>(temp_bodies[i].acceleration) * dt,
            std::get<2>(temp_bodies[i].acceleration) * dt
        );

        std::get<0>(temp_bodies[i].position) = std::get<0>(bodies[i].position) + std::get<0>(k3[i].pos);
        std::get<1>(temp_bodies[i].position) = std::get<1>(bodies[i].position) + std::get<1>(k3[i].pos);
        std::get<2>(temp_bodies[i].position) = std::get<2>(bodies[i].position) + std::get<2>(k3[i].pos);
        std::get<0>(temp_bodies[i].velocity) = std::get<0>(bodies[i].velocity) + std::get<0>(k3[i].vel);
        std::get<1>(temp_bodies[i].velocity) = std::get<1>(bodies[i].velocity) + std::get<1>(k3[i].vel);
        std::get<2>(temp_bodies[i].velocity) = std::get<2>(bodies[i].velocity) + std::get<2>(k3[i].vel);
    }


    calculate_accelerations(temp_bodies);
    for (size_t i = 0; i < bodies.size(); ++i)
    {
        k4[i].pos = std::make_tuple(
            std::get<0>(temp_bodies[i].velocity) * dt,
            std::get<1>(temp_bodies[i].velocity) * dt,
            std::get<2>(temp_bodies[i].velocity) * dt
        );
        k4[i].vel = std::make_tuple(
            std::get<0>(temp_bodies[i].acceleration) * dt,
            std::get<1>(temp_bodies[i].acceleration) * dt,
            std::get<2>(temp_bodies[i].acceleration) * dt
        );

        std::get<0>(bodies[i].position) += (std::get<0>(k1[i].pos) + 3 * std::get<0>(k2[i].pos) + 3 * std::get<0>(k3[i].pos) + std::get<0>(k4[i].pos)) / 8;
        std::get<1>(bodies[i].position) += (std::get<1>(k1[i].pos) + 3 * std::get<1>(k2[i].pos) + 3 * std::get<1>(k3[i].pos) + std::get<1>(k4[i].pos)) / 8;
        std::get<2>(bodies[i].position) += (std::get<2>(k1[i].pos) + 3 * std::get<2>(k2[i].pos) + 3 * std::get<2>(k3[i].pos) + std::get<2>(k4[i].pos)) / 8;

        std::get<0>(bodies[i].velocity) += (std::get<0>(k1[i].vel) + 3 * std::get<0>(k2[i].vel) + 3 * std::get<0>(k3[i].vel) + std::get<0>(k4[i].vel)) / 8;
        std::get<1>(bodies[i].velocity) += (std::get<1>(k1[i].vel) + 3 * std::get<1>(k2[i].vel) + 3 * std::get<1>(k3[i].vel) + std::get<1>(k4[i].vel)) / 8;
        std::get<2>(bodies[i].velocity) += (std::get<2>(k1[i].vel) + 3 * std::get<2>(k2[i].vel) + 3 * std::get<2>(k3[i].vel) + std::get<2>(k4[i].vel)) / 8;
    }
}


void runge_kutta5_step(std::vector<CelestialBody>& bodies, double dt) {
    struct State {
        std::tuple<double, double, double> pos;
        std::tuple<double, double, double> vel;
    };

    const size_t n = bodies.size();
    std::vector<CelestialBody> original = bodies;
    std::vector<State> k1(n), k2(n), k3(n), k4(n), k5(n), k6(n);

    // Этап k1
    calculate_accelerations(original);
    for (size_t i = 0; i < n; ++i) {
        k1[i].pos = original[i].velocity;
        k1[i].vel = original[i].acceleration;
    }

    // Этап k2
    auto temp = original;
    for (size_t i = 0; i < n; ++i) {
        std::get<0>(temp[i].position) += dt * (1.0 / 5.0) * std::get<0>(k1[i].pos);
        std::get<1>(temp[i].position) += dt * (1.0 / 5.0) * std::get<1>(k1[i].pos);
        std::get<2>(temp[i].position) += dt * (1.0 / 5.0) * std::get<2>(k1[i].pos);
        std::get<0>(temp[i].velocity) += dt * (1.0 / 5.0) * std::get<0>(k1[i].vel);
        std::get<1>(temp[i].velocity) += dt * (1.0 / 5.0) * std::get<1>(k1[i].vel);
        std::get<2>(temp[i].velocity) += dt * (1.0 / 5.0) * std::get<2>(k1[i].vel);
    }
    calculate_accelerations(temp);
    for (size_t i = 0; i < n; ++i) {
        k2[i].pos = temp[i].velocity;
        k2[i].vel = temp[i].acceleration;
    }

    // Этап k3
    temp = original;
    for (size_t i = 0; i < n; ++i) {
        std::get<0>(temp[i].position) += dt * (3.0 / 40.0 * std::get<0>(k1[i].pos)) + dt * (9.0 / 40.0 * std::get<0>(k2[i].pos));
        std::get<1>(temp[i].position) += dt * (3.0 / 40.0 * std::get<1>(k1[i].pos)) + dt * (9.0 / 40.0 * std::get<1>(k2[i].pos));
        std::get<2>(temp[i].position) += dt * (3.0 / 40.0 * std::get<2>(k1[i].pos)) + dt * (9.0 / 40.0 * std::get<2>(k2[i].pos));
        std::get<0>(temp[i].velocity) += dt * (3.0 / 40.0 * std::get<0>(k1[i].vel)) + dt * (9.0 / 40.0 * std::get<0>(k2[i].vel));
        std::get<1>(temp[i].velocity) += dt * (3.0 / 40.0 * std::get<1>(k1[i].vel)) + dt * (9.0 / 40.0 * std::get<1>(k2[i].vel));
        std::get<2>(temp[i].velocity) += dt * (3.0 / 40.0 * std::get<2>(k1[i].vel)) + dt * (9.0 / 40.0 * std::get<2>(k2[i].vel));
    }
    calculate_accelerations(temp);
    for (size_t i = 0; i < n; ++i) {
        k3[i].pos = temp[i].velocity;
        k3[i].vel = temp[i].acceleration;
    }

    // Этап k4
    temp = original;
    for (size_t i = 0; i < n; ++i) {
        std::get<0>(temp[i].position) += dt * (44.0 / 45.0 * std::get<0>(k1[i].pos)) - dt * (56.0 / 15.0 * std::get<0>(k2[i].pos)) + dt * (32.0 / 9.0 * std::get<0>(k3[i].pos));
        std::get<1>(temp[i].position) += dt * (44.0 / 45.0 * std::get<1>(k1[i].pos)) - dt * (56.0 / 15.0 * std::get<1>(k2[i].pos)) + dt * (32.0 / 9.0 * std::get<1>(k3[i].pos));
        std::get<2>(temp[i].position) += dt * (44.0 / 45.0 * std::get<2>(k1[i].pos)) - dt * (56.0 / 15.0 * std::get<2>(k2[i].pos)) + dt * (32.0 / 9.0 * std::get<2>(k3[i].pos));
        std::get<0>(temp[i].velocity) += dt * (44.0 / 45.0 * std::get<0>(k1[i].vel)) - dt * (56.0 / 15.0 * std::get<0>(k2[i].vel)) + dt * (32.0 / 9.0 * std::get<0>(k3[i].vel));
        std::get<1>(temp[i].velocity) += dt * (44.0 / 45.0 * std::get<1>(k1[i].vel)) - dt * (56.0 / 15.0 * std::get<1>(k2[i].vel)) + dt * (32.0 / 9.0 * std::get<1>(k3[i].vel));
        std::get<2>(temp[i].velocity) += dt * (44.0 / 45.0 * std::get<2>(k1[i].vel)) - dt * (56.0 / 15.0 * std::get<2>(k2[i].vel)) + dt * (32.0 / 9.0 * std::get<2>(k3[i].vel));
    }
    calculate_accelerations(temp);
    for (size_t i = 0; i < n; ++i) {
        k4[i].pos = temp[i].velocity;
        k4[i].vel = temp[i].acceleration;
    }

    // Этап k5
    temp = original;
    for (size_t i = 0; i < n; ++i) {
        std::get<0>(temp[i].position) += dt * (19372.0 / 6561.0 * std::get<0>(k1[i].pos)) - dt * (25360.0 / 2187.0 * std::get<0>(k2[i].pos)) + dt * (64448.0 / 6561.0 * std::get<0>(k3[i].pos)) - dt * (212.0 / 729.0 * std::get<0>(k4[i].pos));
        std::get<1>(temp[i].position) += dt * (19372.0 / 6561.0 * std::get<1>(k1[i].pos)) - dt * (25360.0 / 2187.0 * std::get<1>(k2[i].pos)) + dt * (64448.0 / 6561.0 * std::get<1>(k3[i].pos)) - dt * (212.0 / 729.0 * std::get<1>(k4[i].pos));
        std::get<2>(temp[i].position) += dt * (19372.0 / 6561.0 * std::get<2>(k1[i].pos)) - dt * (25360.0 / 2187.0 * std::get<2>(k2[i].pos)) + dt * (64448.0 / 6561.0 * std::get<2>(k3[i].pos)) - dt * (212.0 / 729.0 * std::get<2>(k4[i].pos));
        std::get<0>(temp[i].velocity) += dt * (19372.0 / 6561.0 * std::get<0>(k1[i].vel)) - dt * (25360.0 / 2187.0 * std::get<0>(k2[i].vel)) + dt * (64448.0 / 6561.0 * std::get<0>(k3[i].vel)) - dt * (212.0 / 729.0 * std::get<0>(k4[i].vel));
        std::get<1>(temp[i].velocity) += dt * (19372.0 / 6561.0 * std::get<1>(k1[i].vel)) - dt * (25360.0 / 2187.0 * std::get<1>(k2[i].vel)) + dt * (64448.0 / 6561.0 * std::get<1>(k3[i].vel)) - dt * (212.0 / 729.0 * std::get<1>(k4[i].vel));
        std::get<2>(temp[i].velocity) += dt * (19372.0 / 6561.0 * std::get<2>(k1[i].vel)) - dt * (25360.0 / 2187.0 * std::get<2>(k2[i].vel)) + dt * (64448.0 / 6561.0 * std::get<2>(k3[i].vel)) - dt * (212.0 / 729.0 * std::get<2>(k4[i].vel));
    }
    calculate_accelerations(temp);
    for (size_t i = 0; i < n; ++i) {
        k5[i].pos = temp[i].velocity;
        k5[i].vel = temp[i].acceleration;
    }

    // Этап k6
    temp = original;
    for (size_t i = 0; i < n; ++i) {
        std::get<0>(temp[i].position) += dt * (9017.0 / 3168.0 * std::get<0>(k1[i].pos)) - dt * (355.0 / 33.0 * std::get<0>(k2[i].pos)) + dt * (46732.0 / 5247.0 * std::get<0>(k3[i].pos)) + dt * (49.0 / 176.0 * std::get<0>(k4[i].pos)) - dt * (5103.0 / 18656.0 * std::get<0>(k5[i].pos));
        std::get<1>(temp[i].position) += dt * (9017.0 / 3168.0 * std::get<1>(k1[i].pos)) - dt * (355.0 / 33.0 * std::get<1>(k2[i].pos)) + dt * (46732.0 / 5247.0 * std::get<1>(k3[i].pos)) + dt * (49.0 / 176.0 * std::get<1>(k4[i].pos)) - dt * (5103.0 / 18656.0 * std::get<1>(k5[i].pos));
        std::get<2>(temp[i].position) += dt * (9017.0 / 3168.0 * std::get<2>(k1[i].pos)) - dt * (355.0 / 33.0 * std::get<2>(k2[i].pos)) + dt * (46732.0 / 5247.0 * std::get<2>(k3[i].pos)) + dt * (49.0 / 176.0 * std::get<2>(k4[i].pos)) - dt * (5103.0 / 18656.0 * std::get<2>(k5[i].pos));
        std::get<0>(temp[i].velocity) += dt * (9017.0 / 3168.0 * std::get<0>(k1[i].vel)) - dt * (355.0 / 33.0 * std::get<0>(k2[i].vel)) + dt * (46732.0 / 5247.0 * std::get<0>(k3[i].vel)) + dt * (49.0 / 176.0 * std::get<0>(k4[i].vel)) - dt * (5103.0 / 18656.0 * std::get<0>(k5[i].vel));
        std::get<1>(temp[i].velocity) += dt * (9017.0 / 3168.0 * std::get<1>(k1[i].vel)) - dt * (355.0 / 33.0 * std::get<1>(k2[i].vel)) + dt * (46732.0 / 5247.0 * std::get<1>(k3[i].vel)) + dt * (49.0 / 176.0 * std::get<1>(k4[i].vel)) - dt * (5103.0 / 18656.0 * std::get<1>(k5[i].vel));
        std::get<2>(temp[i].velocity) += dt * (9017.0 / 3168.0 * std::get<2>(k1[i].vel)) - dt * (355.0 / 33.0 * std::get<2>(k2[i].vel)) + dt * (46732.0 / 5247.0 * std::get<2>(k3[i].vel)) + dt * (49.0 / 176.0 * std::get<2>(k4[i].vel)) - dt * (5103.0 / 18656.0 * std::get<2>(k5[i].vel));
    }
    calculate_accelerations(temp);
    for (size_t i = 0; i < n; ++i) {
        k6[i].pos = temp[i].velocity;
        k6[i].vel = temp[i].acceleration;
    }

    // Финальное обновление
    for (size_t i = 0; i < n; ++i) {
        const double b1 = 35.0 / 384.0, b3 = 500.0 / 1113.0, b4 = 125.0 / 192.0, b5 = -2187.0 / 6784.0, b6 = 11.0 / 84.0;
        std::get<0>(bodies[i].position) = std::get<0>(original[i].position) + dt * (b1 * std::get<0>(k1[i].pos) + b3 * std::get<0>(k3[i].pos) + b4 * std::get<0>(k4[i].pos) + b5 * std::get<0>(k5[i].pos) + b6 * std::get<0>(k6[i].pos));
        std::get<1>(bodies[i].position) = std::get<1>(original[i].position) + dt * (b1 * std::get<1>(k1[i].pos) + b3 * std::get<1>(k3[i].pos) + b4 * std::get<1>(k4[i].pos) + b5 * std::get<1>(k5[i].pos) + b6 * std::get<1>(k6[i].pos));
        std::get<2>(bodies[i].position) = std::get<2>(original[i].position) + dt * (b1 * std::get<2>(k1[i].pos) + b3 * std::get<2>(k3[i].pos) + b4 * std::get<2>(k4[i].pos) + b5 * std::get<2>(k5[i].pos) + b6 * std::get<2>(k6[i].pos));

        std::get<0>(bodies[i].velocity) = std::get<0>(original[i].velocity) + dt * (b1 * std::get<0>(k1[i].vel) + b3 * std::get<0>(k3[i].vel) + b4 * std::get<0>(k4[i].vel) + b5 * std::get<0>(k5[i].vel) + b6 * std::get<0>(k6[i].vel));
        std::get<1>(bodies[i].velocity) = std::get<1>(original[i].velocity) + dt * (b1 * std::get<1>(k1[i].vel) + b3 * std::get<1>(k3[i].vel) + b4 * std::get<1>(k4[i].vel) + b5 * std::get<1>(k5[i].vel) + b6 * std::get<1>(k6[i].vel));
        std::get<2>(bodies[i].velocity) = std::get<2>(original[i].velocity) + dt * (b1 * std::get<2>(k1[i].vel) + b3 * std::get<2>(k3[i].vel) + b4 * std::get<2>(k4[i].vel) + b5 * std::get<2>(k5[i].vel) + b6 * std::get<2>(k6[i].vel));
    }
}


void save_trajectories(std::vector<TrajectoryPoint> trajectories, std::string filename)
{
    std::ofstream file(filename);
    file << "Time(h),Body,X(m),Y(m),Z(m)\n";
    for(auto& point : trajectories)
    {
        file << point.time << ","
            << point.name << ","
            << std::get<0>(point.pos) << ","
            << std::get<1>(point.pos) << ","
            << std::get<2>(point.pos) << "\n";
    }
    file.close();
}

int main()
{
    std::vector<CelestialBody> bodies = {
        // Солнце
        {"Sun", -1.9885e30, {0, 0, 0}, {0, 0, 0}, {0, 0, 0}},

        // Меркурий
        {"Mercury", 3.3011e23, {57.9e9, 0, 0}, {0, 47.36e3, 0}, {0, 0, 0}},

        // Венера
        {"Venus", 4.8675e24, {108.2e9, 0, 0}, {0, 35.02e3, 0}, {0, 0, 0}},

        // Земля
        {"Earth", 5.9736e24, {1.496e11, 0, 0}, {0, 29.78e3, 0}, {0, 0, 0}},

        // Луна
        {"Moon", 7.342e22, {1.496e11 + 3.844e8, 0, 0}, {0, 29.78e3 + 1.022e3, 0}, {0, 0, 0}},

        // Марс
        {"Mars", 6.4171e23, {227.9e9, 0, 0}, {0, 24.07e3, 0}, {0, 0, 0}},

        // Юпитер
        {"Jupiter", 1.8982e27, {778.5e9, 0, 0}, {0, 13.07e3, 0}, {0, 0, 0}},

        // Сатурн
        {"Saturn", 5.6834e26, {1.433e12, 0, 0}, {0, 9.69e3, 0}, {0, 0, 0}},

        // Уран
        {"Uranus", 8.6810e25, {2.872e12, 0, 0}, {0, 6.81e3, 0}, {0, 0, 0}},
       
        // Нептун
        {"Neptune", 1.0241e26, {4.495e12, 0, 0}, {0, 5.43e3, 0}, {0, 0, 0}},

        // Плутон
        {"Pluto", 1.303e22, {5.906e12, 0, 0}, {0, 4.67e3, 0}, {0, 0, 0}}
    };

    std::vector<TrajectoryPoint> trajectories;
    double dt = 3600;

    int choice = 0;
    std::cout << "option 1 = basic calculations" << std::endl;
    std::cout << "option 2 = no sun" << std::endl;
    std::cout << "option 3 = no earth" << std::endl;
    std::cout << "option 4 = sun has speed" << std::endl;
    std::cout << "input option (1/2/3/4): ";
    std::cin >> choice;

    switch (choice)
    {
    case 2:
        bodies[0].mass = 1;
        break;

    case 3:
        bodies[3].mass = 1;
        break;

    case 4:

        bodies[0].velocity = { 13.07e3, 0, 13.07e3 };
        break;

    case 1:
        break;

    default:
        std::cout << "wtf?" << std::endl;
        return 0;
    }

    for (double gt = 0; gt <= 1000000 * 24 * 3600; gt += 3600)
    {
        runge_kutta_step(bodies, dt);

        for (auto body : bodies)
        {
            trajectories.push_back({ gt, body.name, body.position });
        }
    }

    save_trajectories(trajectories, "trajectories3.csv");
    std::cout << "vaaaaaaaaaagh!" << std::endl;
  

    return 0;
}