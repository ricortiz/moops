#ifndef WRITE_VTK
#define WRITE_VTK
/// C++ Interfaces for writing vtk
///
/// @brief Description: Write vtk xml formated data
///
/// Author: Ricardo Ortiz <ortiz@unc.edu>, (C) 2008
/// $Id $
#include <iostream>
#include <sstream>
#include <cassert>
#include <vtkXMLUnstructuredGridWriter.h>
#include <vtkUnstructuredGrid.h>
#include <vtkPointData.h>
#include <vtkCellArray.h>
#include <vtkDoubleArray.h>
#include <vtkSmartPointer.h>
#include <vtkXMLPolyDataWriter.h>
#include <vtksys/SystemTools.hxx>

#include<cmath>

namespace IO
{
    namespace detail
    {
        inline std::string file_number(size_t counter, size_t total_digits = 6)
        {
            std::stringstream filename;

            size_t digits = (size_t)std::floor(counter > 9 ? std::log10(counter) + 1 : 1);
            assert(digits <= total_digits || !"Number of digits is bigger than allowed.");
            std::string zeros(total_digits - digits, '0');

            // concatenate zeros + counter
            filename << zeros << counter;
            return filename.str();
        }
    }
    template<typename vtk_writer_type>
    class VtkWriter
    {
            size_t m_file_counter;
            std::string m_data_path;
            vtkSmartPointer<vtk_writer_type> m_writer;
            bool m_write_binary;

        public:
            VtkWriter() : m_file_counter(0), m_data_path("./"), m_writer(vtkSmartPointer<vtk_writer_type>::New()), m_write_binary(false)  {}
            VtkWriter(const std::string &data_path, bool write_binary = false) :
                    m_file_counter(0),
                    m_data_path(data_path),
                    m_writer(vtkSmartPointer<vtk_writer_type>::New()),
                    m_write_binary(write_binary)
            {
                if (m_write_binary)
                    m_writer->SetDataModeToBinary();
                else
                    m_writer->SetDataModeToAscii();
            }

            template<typename input_type>
            void write(std::vector<input_type> &poly_data)
            {
                for(size_t i = 0; i < poly_data.size(); ++i)
                    write(poly_data[i],0);
            }
            
            template<typename data_set_type>
            void write(data_set_type &data_set, double timestep, bool set_number = true, bool print = true)
            {
                m_writer->SetInput(data_set);
                std::string file_name = m_data_path;
                if(set_number) file_name += "_" + detail::file_number(m_file_counter++, 10);
                file_name += ".";
                file_name += m_writer->GetDefaultFileExtension();
                m_writer->SetFileName(file_name.c_str());
                if (print) std::cout << "Saving " << file_name << " ... ";
                m_writer->WriteNextTime(timestep);
                if (print) std::cout << "done." << std::endl;
            }

            void setDataPath(std::string data_path)
            {
                m_data_path = data_path;
            }

    };

}
#endif




