#ifndef WRITE_VTK
#define WRITE_VTK
/// C++ Interfaces for writing vtk
///
/// @brief Description: Write vtk xml formated data
///
/// Author: Ricardo Ortiz <ortiz@unc.edu>, (C) 2008
/// $Id $
#include <iostream>
#include <vtkXMLUnstructuredGridWriter.h>
#include <vtkUnstructuredGrid.h>
#include <vtkPointData.h>
#include <vtkCellArray.h>
#include <vtkDoubleArray.h>
#include <vtkSmartPointer.h>

#include <vtksys/SystemTools.hxx>

namespace IO
{

    template<typename vtk_writer_type>
    class VTKWriter
    {
            size_t m_file_counter;
            std::string m_data_path;
            vtkSmartPointer<vtk_writer_type> m_writer;

        public:
            VTKWriter() : m_file_counter(0){}
            VTKWriter(const std::string &data_path) :
                    m_file_counter(0),
                    m_data_path(data_path),
                    m_writer(vtkSmartPointer<vtk_writer_type>::New())
            {
//                 m_writer->SetDataModeToAscii();
                m_writer->SetDataModeToBinary();
            }    

            template<typename input_type>
            void setInput(input_type &data)
            {
                m_writer->SetInput(data);
            }
            
            void write(double timestep, bool print = true)
            {
                std::string file_name = m_data_path + "data_" + file_number(m_file_counter++,10) + "." + m_writer->GetDefaultFileExtension();
                m_writer->SetFileName(file_name.c_str());
                if(print) std::cout << "Saving " << file_name << " ... ";
                m_writer->WriteNextTime(timestep);
                if(print) std::cout << "done." << std::endl;
            }

            inline std::string file_number(size_t counter, int total_digits = 6)
            {
                std::stringstream filename;

                size_t digits = (size_t)std::floor(counter > 9 ? std::log10(counter) + 1 : 1);
                assert(digits <= total_digits || !"Number of digits is bigger than allowed.");
                std::string zeros(total_digits - digits,'0');

                // concatenate zeros + counter
                filename << zeros << counter;
                return filename.str();
            }
            
    };
}
#endif




