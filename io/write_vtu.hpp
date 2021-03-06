#ifndef WRITE_VTK
#define WRITE_VTK
/****************************************************************************
** MOOPS -- Modular Object Oriented Particle Simulator
** Copyright (C) 2011-2012  Ricardo Ortiz <ortiz@unc.edu>
**
** This program is free software: you can redistribute it and/or modify
** it under the terms of the GNU General Public License as published by
** the Free Software Foundation, either version 3 of the License, or
** (at your option) any later version.
**
** This program is distributed in the hope that it will be useful,
** but WITHOUT ANY WARRANTY; without even the implied warranty of
** MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
** GNU General Public License for more details.
**
** You should have received a copy of the GNU General Public License
** along with this program.  If not, see <http://www.gnu.org/licenses/>.
****************************************************************************/
#include <iostream>
#include <sstream>
#include <cmath>
#include <stdexcept>

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
            if(digits > total_digits)
                throw std::range_error("file_number(): Number of digits is bigger than allowed.");
            std::string zeros(total_digits - digits, '0');

            // concatenate zeros + counter
            filename << zeros << counter;
            return filename.str();
        }
    }

    /** \brief Wrapper for VTK generic writer
     *
     *  \tparam vtk_writer_type is the VTK's writer
     *
     */
    template<typename vtk_writer_type>
    class VtkWriter
    {
            size_t m_file_counter;
            std::string m_data_path;
            vtkSmartPointer<vtk_writer_type> m_writer;
            bool m_write_binary;

        public:
            VtkWriter() : m_file_counter(0), m_data_path("./"), m_writer(vtkSmartPointer<vtk_writer_type>::New()), m_write_binary(true)  {}
            VtkWriter(const std::string &data_path, bool write_binary = true) :
                    m_file_counter(0),
                    m_data_path(data_path),
                    m_writer(vtkSmartPointer<vtk_writer_type>::New()),
                    m_write_binary(write_binary)
            {
                vtksys::SystemTools::MakeDirectory(data_path.c_str());
                if (m_write_binary)
                    m_writer->SetDataModeToBinary();
                else
                    m_writer->SetDataModeToAscii();
            }
            
            /**
             * \brief Convenience method for writing many datasets into files
             *
             * \param poly_data vector of datasets
             * 
             **/
            template<typename input_type>
            void write(std::vector<input_type> &data)
            {
                for(size_t i = 0; i < data.size(); ++i)
                    write(data[i],0,1);
            }
            
            /**
             * \brief Writes the dataset into a file
             *
             * \param data_set data
             * \param time current time
             *
             **/
            template<typename data_set_type>
            void write(data_set_type &data_set, double time, bool set_number = true, bool print = true)
            {
                m_writer->SetInput(data_set);
                std::string file_name = m_data_path;
                if(set_number) file_name += "_" + detail::file_number(m_file_counter++, 10);
                file_name += ".";
                file_name += m_writer->GetDefaultFileExtension();
                m_writer->SetFileName(file_name.c_str());
                if (print) std::cout << "Saving " << file_name << " ... ";
                m_writer->WriteNextTime(time);
                if (print) std::cout << "done." << std::endl;
            }

            /**
             * \brief Sets the datapath where the files are going to be saved
             *
             * \param data_path string containing the path to the data directory
             *
             **/
            void setDataPath(const std::string &data_path)
            {
                vtksys::SystemTools::MakeDirectory(data_path.c_str());
                m_data_path = data_path;
            }

    };

}
#endif




