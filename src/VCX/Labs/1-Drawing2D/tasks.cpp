#include <cstdlib>
#include <random>

#include <spdlog/spdlog.h>

#include "Labs/1-Drawing2D/tasks.h"

using VCX::Labs::Common::ImageRGB;

namespace VCX::Labs::Drawing2D {
    /******************* 1.Image Dithering *****************/
    void DitheringThreshold(
        ImageRGB &       output,
        ImageRGB const & input) {
        for (std::size_t x = 0; x < input.GetSizeX(); ++x)
            for (std::size_t y = 0; y < input.GetSizeY(); ++y) {
                glm::vec3 color = input.At(x, y);
                output.At(x, y) = {
                    color.r > 0.5 ? 1 : 0,
                    color.g > 0.5 ? 1 : 0,
                    color.b > 0.5 ? 1 : 0,
                };
            }
    }

    void DitheringRandomUniform(
        ImageRGB &       output,
        ImageRGB const & input) {
        std::random_device rd;
        std::mt19937       gen(rd()); // Mersenne Twister算法
        std::uniform_real_distribution<float> dis(-0.5, 0.5);
        for(std::size_t x = 0; x < input.GetSizeX(); ++x)
            for(std::size_t y = 0; y < input.GetSizeY(); ++y) {
                glm::vec3 color = input.At(x, y);

                float perturb = dis(gen);


                output.At(x, y) = {
                    color.r+perturb > 0.5 ? 1 : 0,
                    color.g+perturb > 0.5 ? 1 : 0,
                    color.b+perturb > 0.5 ? 1 : 0,
                };
            }
        // your code here:
    }

    void DitheringRandomBlueNoise(
        ImageRGB &       output,
        ImageRGB const & input,
        ImageRGB const & noise) {
        
        for (std::size_t x = 0; x < input.GetSizeX(); ++x)
            for (std::size_t y = 0; y < input.GetSizeY(); ++y) {
                glm::vec3 color = input.At(x, y);
                glm::vec3 noiseColor = noise.At(x, y);
                output.At(x, y) = {
                    color.r+noiseColor.r >1  ? 1 : 0,
                    color.g+noiseColor.g >1  ? 1 : 0,
                    color.b+noiseColor.b >1  ? 1 : 0,
                };
            }
        // your code here:
    }

    void DitheringOrdered(
        ImageRGB &       output,
        ImageRGB const & input) {
        const int dithering_matirx[3][3]={
            {6, 8, 4},
            {1, 0, 3},
            {5, 2, 7}
        };
        // output.resize(input.GetSizeX()*3,input.GetSizeY()*3);
        for (int y=0;y<input.GetSizeY();y++){
            for (int x=0;x<input.GetSizeX();x++){
                glm::vec3 color = input.At(x, y);
                for (int i=0;i<3;i++){
                    for (int j=0;j<3;j++){
                        output.At(x*3+i, y*3+j) = {
                            color.r > (dithering_matirx[i][j]/9.0) ? 1 : 0,
                            color.g > (dithering_matirx[i][j]/9.0) ? 1 : 0,
                            color.b > (dithering_matirx[i][j]/9.0) ? 1 : 0,
                        };
                    }
                }
            }
        }
        // your code here:
    }

    void DitheringErrorDiffuse(
        ImageRGB &       output,
        ImageRGB const & input) {
        ImageRGB tmp_input=input;
        for (int y=0;y<input.GetSizeY();y++){
            for (int x=0;x<input.GetSizeX();x++){
                glm::vec3 color = tmp_input.At(x, y);
               output.At(x, y) = {
                    color.r >0.5? 1 : 0,
                    color.g >0.5? 1 : 0,
                    color.b >0.5? 1 : 0,
               };
               glm::vec3 the_color = output.At(x, y);
               double error=color.r-the_color.r;
               if(x<input.GetSizeX()-1){
                glm::vec3 tmp_color = tmp_input.At(x+1, y);
                tmp_input.At(x+1, y) = {
                    tmp_color.r+error*0.4375,
                    tmp_color.g+error*0.4375,
                    tmp_color.b+error*0.4375,
                };
               };
               if(y<input.GetSizeY()-1){
                    if (x>0){
                        glm::vec3 tmp_color = tmp_input.At(x-1, y+1);
                        tmp_input.At(x-1, y+1) = {
                            tmp_color.r+error*0.1875,
                            tmp_color.g+error*0.1875,
                            tmp_color.b+error*0.1875,
                        };
                    }
                    if (x<input.GetSizeX()-1){
                        glm::vec3 tmp_color = tmp_input.At(x+1, y+1);
                        tmp_input.At(x+1, y+1) = {
                            tmp_color.r+error*0.0625,
                            tmp_color.g+error*0.0625,
                            tmp_color.b+error*0.0625,
                        };
                    }
                glm::vec3 tmp_color = tmp_input.At(x, y+1);
                tmp_input.At(x, y+1) = {
                    tmp_color.r+error*0.3125,
                    tmp_color.g+error*0.3125,
                    tmp_color.b+error*0.3125,
                };
               };


            };
        };
        // your code here:
    }

    /******************* 2.Image Filtering *****************/
    void Blur(
        ImageRGB &       output,
        ImageRGB const & input) {
        const int kernel[3][3] = {
            {1, 1, 1},
            {1, 1, 1},
            {1, 1, 1},
        };
        for (int y = 0; y < input.GetSizeY(); y++) {
            for (int x = 0; x < input.GetSizeX(); x++) {
                glm::vec3 color = {0, 0, 0};
                for (int i = 0; i < 3; i++) {
                    for (int j = 0; j < 3; j++) {
                        int new_x = x + i - 1;
                        int new_y = y + j - 1;
                        if (new_x >= 0 && new_x < input.GetSizeX() && new_y >= 0 && new_y < input.GetSizeY()) {
                            glm::vec3 tmp_color = input.At(new_x, new_y);
                            color.r+=tmp_color.r*kernel[i][j]/9.0;
                            color.g+=tmp_color.g*kernel[i][j]/9.0;
                            color.b+=tmp_color.b*kernel[i][j]/9.0;
                        }
                    }
                }
                output.At(x, y) = color;
            }
        }
        // your code here:
       
    }

    void Edge(
        ImageRGB &       output,
        ImageRGB const & input) {
        const int kernel[3][3] = {
            {1, 1, 1},
            {1, -8, 1},
            {1, 1, 1},
        };
        for (int y = 0; y < input.GetSizeY(); y++) {
            for (int x = 0; x < input.GetSizeX(); x++) {
                glm::vec3 color = {0, 0, 0};
                for (int i = 0; i < 3; i++) {
                    for (int j = 0; j < 3; j++) {
                        int new_x = x + i - 1;
                        int new_y = y + j - 1;
                        if (new_x >= 0 && new_x < input.GetSizeX() && new_y >= 0 && new_y < input.GetSizeY()) {
                            glm::vec3 tmp_color = input.At(new_x, new_y);
                            color.r+=tmp_color.r*kernel[i][j];
                            color.g+=tmp_color.g*kernel[i][j];
                            color.b+=tmp_color.b*kernel[i][j];
                        }
                    }
                }
                output.At(x, y) = color;
            }
        }
        // your code here:
    }

    /******************* 3. Image Inpainting *****************/
    void Inpainting(
        ImageRGB &         output,
        ImageRGB const &   inputBack,
        ImageRGB const &   inputFront,
        const glm::ivec2 & offset) {
        output             = inputBack;
        std::size_t width  = inputFront.GetSizeX();
        std::size_t height = inputFront.GetSizeY();
        glm::vec3 * g      = new glm::vec3[width * height];
        memset(g, 0, sizeof(glm::vec3) * width * height);
        // set boundary condition
        for (std::size_t y = 0; y < height; ++y) {
            // set boundary for (0, y), your code: g[y * width] = ?
            g[y * width] =  inputBack.At(offset.x, y + offset.y)-inputFront.At(0, y);
            // set boundary for (width - 1, y), your code: g[y * width + width - 1] = ?
            g[y * width + width - 1] =  inputBack.At(offset.x + width - 1, y + offset.y)-inputFront.At(width - 1, y);
        }
        for (std::size_t x = 0; x < width; ++x) {
            // set boundary for (x, 0), your code: g[x] = ?
            g[x] =  inputBack.At(x + offset.x, offset.y)-inputFront.At(x, 0);
            // set boundary for (x, height - 1), your code: g[(height - 1) * width + x] = ?
            g[(height - 1) * width + x] =  inputBack.At(x + offset.x, offset.y + height - 1)-inputFront.At(x, height - 1);
        }

        // Jacobi iteration, solve Ag = b
        for (int iter = 0; iter < 8000; ++iter) {
            for (std::size_t y = 1; y < height - 1; ++y)
                for (std::size_t x = 1; x < width - 1; ++x) {
                    g[y * width + x] = (g[(y - 1) * width + x] + g[(y + 1) * width + x] + g[y * width + x - 1] + g[y * width + x + 1]);
                    g[y * width + x] = g[y * width + x] * glm::vec3(0.25);
                }
        }

        for (std::size_t y = 0; y < inputFront.GetSizeY(); ++y)
            for (std::size_t x = 0; x < inputFront.GetSizeX(); ++x) {
                glm::vec3 color = g[y * width + x] + inputFront.At(x, y);
                output.At(x + offset.x, y + offset.y) = color;
                
            }
        delete[] g;
    }

    /******************* 4. Line Drawing *****************/
    void DrawLine(
        ImageRGB &       canvas,
        glm::vec3 const  color,
        glm::ivec2 const p0,
        glm::ivec2 const p1) {
        int x0=p0.x,y0=p0.y,x1=p1.x,y1=p1.y;
        bool steep=std::abs(y1-y0)>std::abs(x1-x0);
        if (steep){
            std::swap(x0,y0);
            std::swap(x1,y1);
        }
        if (x0>x1){
            std::swap(x0,x1);
            std::swap(y0,y1);
        }
        int dx=x1-x0,dy=std::abs(y1-y0);
        int error=dx/2;
        int ystep=(y0<y1)?1:-1;
        int y=y0;
        for(int x=x0;x<=x1;x++){
            if (steep){
                canvas.At(y,x)=color;
            }
            else{
                canvas.At(x,y)=color;
            }
            error-=dy;
            if(error<0){
                y+=ystep;
                error+=dx;
            }
        }
        // your code here:
    }

    /******************* 5. Triangle Drawing *****************/
    void DrawTriangleFilled(
        ImageRGB &       canvas,
        glm::vec3 const  color,
        glm::ivec2 const p0,
        glm::ivec2 const p1,
        glm::ivec2 const p2) {
        int x0=p0.x,y0=p0.y,x1=p1.x,y1=p1.y,x2=p2.x,y2=p2.y;
        if(y1<y0){
            std::swap(x0,x1);
            std::swap(y0,y1);
        }
        if(y2<y0){
            std::swap(x0,x2);
            std::swap(y0,y2);
        }
        if(y2<y1){
            std::swap(x1,x2);
            std::swap(y1,y2);
        }
        //y0<=y1<=y2
        
        auto interpolate=[](int x1,int y1,int x2,int y2,int y){
            if(y1==y2){
                return x1;
            }
            return x1+(x2-x1)*(y-y1)/(y2-y1);
        };

        for(int y=y0;y<=y2;y++){
            bool lowerHalf = y <= y1;
            int x_start=interpolate(x0,y0,x2,y2,y);
            int x_end=lowerHalf?interpolate(x0,y0,x1,y1,y):interpolate(x1,y1,x2,y2,y);
            if(x_start>x_end){
                std::swap(x_start,x_end);
            }
            for(int x=x_start;x<=x_end;x++){
                canvas.At(x,y)=color;
            }
        }
        // your code here:
    }

    /******************* 6. Image Supersampling *****************/
    void Supersample(
        ImageRGB &       output,
        ImageRGB const & input,
        int              rate) {
        for (int y = 0; y < output.GetSizeY(); y++) {
            for (int x = 0; x < output.GetSizeX(); x++) {
                glm::vec3 color = {0, 0, 0};
                for (int i = 0; i < rate; i++) {
                    for (int j = 0; j < rate; j++) {
                        int new_x = x * rate + i;
                        int new_y = y * rate + j;
                        if(new_x<input.GetSizeX()&&new_y<input.GetSizeY()){
                            glm::vec3 tmp_color = input.At(new_x, new_y);
                            color.r+=tmp_color.r;
                            color.g+=tmp_color.g;
                            color.b+=tmp_color.b;
                    }
                }
                output.At(x, y) ={
                    color.r/(rate * rate),
                    color.g/(rate * rate),
                    color.b/(rate * rate),
                    };
                }
            }
        // your code here:
        }
    }
    /******************* 7. Bezier Curve *****************/
    // Note: Please finish the function [DrawLine] before trying this part.
    glm::vec2 CalculateBezierPoint(
        std::span<glm::vec2> points,
        float const          t) {
        if(points.size()==1){
            return points[0];
        }
        std::vector<glm::vec2> new_points;
        for(int i=0;i<points.size()-1;i++){
            new_points.push_back(points[i]*(1-t)+points[i+1]*t);

        }
        return CalculateBezierPoint(new_points,t);
        // your code here:
        return glm::vec2 {0, 0};
    }
} // namespace VCX::Labs::Drawing2D