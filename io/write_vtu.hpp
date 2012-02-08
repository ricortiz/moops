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
#include <vtkXMLPolyDataWriter.h>
#include <vtksys/SystemTools.hxx>

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
    template<typename vtk_storage, typename vtk_writer_type>
    class VtkWriter;

    template<typename vtk_storage>
    class VtkWriter<vtk_storage, vtkXMLPolyDataWriter>
    {
            size_t m_file_counter;
            std::string m_data_path;
            vtkSmartPointer<vtkXMLPolyDataWriter> m_writer;
            bool m_write_binary;

        public:
            VtkWriter() : m_file_counter(0), m_write_binary(true) {}
            VtkWriter(const std::string &data_path, vtk_storage &storage, bool write_binary = true) :
                    m_file_counter(0),
                    m_data_path(data_path),
                    m_writer(vtkSmartPointer<vtkXMLPolyDataWriter>::New()),
                    m_write_binary(write_binary)
            {
                m_writer->SetInput(storage.grid());
                if (m_write_binary)
                    m_writer->SetDataModeToBinary();
                else
                    m_writer->SetDataModeToAscii();
            }

            template<typename input_type>
            void setInput(input_type &poly_data, int i = 0)
            {
                m_writer->SetInput(poly_data, i);
            }

            void write(double timestep, bool print = true)
            {
                std::string file_name = m_data_path + "data_" + detail::file_number(m_file_counter++, 10) + "." + m_writer->GetDefaultFileExtension();
                m_writer->SetFileName(file_name.c_str());
                if (print) std::cout << "Saving " << file_name << " ... ";
                m_writer->WriteNextTime(timestep);
                if (print) std::cout << "done." << std::endl;
            }

    };

}
#endif




