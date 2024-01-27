#define CL_HPP_ENABLE_EXCEPTIONS

#include <CL/opencl.hpp>
#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <random>
#include <SFML/Graphics.hpp>
#include "imgui-master/imgui.h"
#include "imgui-sfml-2.6.x/imgui-SFML.h"

using namespace std;

string load_kernel_file(){
    string filePath = "../kernel.cl";
    ifstream file(filePath);

    if (!file.is_open()) {
        std::cerr << "Error opening file: " << filePath << std::endl;
        exit(1);
    }

    std::string fileContent((std::istreambuf_iterator<char>(file)), std::istreambuf_iterator<char>());
    file.close();

    return fileContent;
}

cl::Device find_device(){
    cl::Device device;

    std::vector<cl::Platform> platforms;
    cl::Platform::get(&platforms);

    for (const auto& p : platforms) {
        std::vector<cl::Device> devices;
        p.getDevices(CL_DEVICE_TYPE_GPU, &devices);

        for (const auto& d : devices) {
            if(d.getInfo<CL_DEVICE_NAME>()=="gfx1030"){
                return d;
            }
        }
    }
    std::cerr << "found no gpu" << std::endl;
    exit(1);
}

std::vector<unsigned char> convertFloatsToBytes(const std::vector<float>& floatValues) {
    std::vector<unsigned char> byteValues;
    byteValues.reserve(floatValues.size());
    std::transform(floatValues.begin(), floatValues.end(), std::back_inserter(byteValues),[](float floatValue) {
                       return static_cast<unsigned char>(std::min(std::max(floatValue * 255.0f, 0.0f), 255.0f));
                   });
    return byteValues;

}

const int width = 2048;
const int height = 1360;

void compute_ray(cl::Kernel& kernel,cl::CommandQueue& queue,cl::Image2D& resultBuffer,std::vector<float>& image_data,double r0,double theta0,double a){
    try {
        cl::Event kernelEvent;
        std::array<cl::size_type, 3> origin{0, 0, 0};
        std::array<cl::size_type, 3> region{width, height, 1};

        //seed
        std::random_device rd;
        std::mt19937 mt(rd());
        std::uniform_int_distribution<int> dist(0, 1000000);
        int random_number = dist(mt);

        kernel.setArg(0, resultBuffer);
        kernel.setArg(1, r0);
        kernel.setArg(2, theta0);
        kernel.setArg(3, a);
        kernel.setArg(4,random_number);

        cl::NDRange globalSize(width, height);
        cl::NDRange localSize(16, 16);

        queue.enqueueNDRangeKernel(kernel, cl::NullRange, globalSize, localSize, nullptr, &kernelEvent);
        kernelEvent.wait();

        queue.enqueueReadImage(resultBuffer, CL_TRUE, origin, region, 0, 0, image_data.data(), nullptr, &kernelEvent);
        kernelEvent.wait();
    }catch (const cl::Error& err) {
        std::cerr << "OpenCL error: " << err.what() << " (" << err.err() << ")" << std::endl;
        return;
    }
}

int main() {


    cl::Device device = find_device();

    cl::Context context(device);
    cl::CommandQueue queue(context);

    string source = load_kernel_file();

    cl::Program program(context, source);

    cl::ImageFormat format(CL_RGBA, CL_FLOAT);
    cl::Image2D resultBuffer(context, CL_MEM_WRITE_ONLY, format, width, height);

    std::vector<float> image_data(width * height * 4);

    cl::Event kernelEvent;
    cl::Kernel compute_ray_kernel;
    try {
        program.build();
        compute_ray_kernel = cl::Kernel(program, "compute_ray");
    }catch (const cl::Error& err) {
        std::cerr << "OpenCL error: " << err.what() << " (" << err.err() << ")" << std::endl;
        std::cerr << "Build log: " << program.getBuildInfo<CL_PROGRAM_BUILD_LOG>(device) << std::endl;
        return 1;
    }



    // Create an SFML window
    sf::RenderWindow window(sf::VideoMode(width, height), "SFML Image Display",sf::Style::Default, sf::ContextSettings(0, 0, 8));

    ImGui::SFML::Init(window);

    sf::Sprite sprite;
    sf::Image image;
    sf::Texture texture;
    while (window.isOpen()) {
        sf::Event event;
        while (window.pollEvent(event)) {
            ImGui::SFML::ProcessEvent(window,event);
            if (event.type == sf::Event::Closed) {
                window.close();
            }
        }

        ImGui::SFML::Update(window, sf::seconds(1.f / 60.f));
        ImGui::Begin("Settings");

        static float r0 = 16.0f;
        static float theta0 = 86.0f;
        static float a = 0.0f;
        ImGui::SliderFloat("Radius", &r0, 0.0f, 100.0f);
        ImGui::SliderFloat("angle", &theta0, 0.01f, 89.9f);
        ImGui::SliderFloat("angular momentum", &a, 0.0f, 1.0f);

        if (ImGui::Button("render")) {
            compute_ray(compute_ray_kernel,queue,resultBuffer,image_data,(double)r0,3.1415*((double)theta0)/180.0,(double)a);
            std::vector<unsigned char> byte_data = convertFloatsToBytes(image_data);
            image.create(width, height, byte_data.data());
            texture.loadFromImage(image);
            sprite = sf::Sprite(texture);
        }

        ImGui::End();
        window.clear();
        window.draw(sprite);
        ImGui::SFML::Render(window);
        window.display();
    }
    ImGui::SFML::Shutdown();

    return 0;

}