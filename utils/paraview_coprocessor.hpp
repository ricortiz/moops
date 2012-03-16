#ifndef PARAVIEW_COPROCESSOR_HPP
#define PARAVIEW_COPROCESSOR_HPP
//=========================================================================
//
//  Program:   Modular Object Oriented Particle Simulator
//  Module:    vtkParticleSystemStorage
//
//  Copyright (c) Ricardo Ortiz
//  All rights reserved.
//     This software is distributed WITHOUT ANY WARRANTY; without even
//     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
//     PURPOSE.
//
//=========================================================================

#include<vtkSmartPointer.h>
#include<vtkCPDataDescription.h>
#include<vtkCPInputDataDescription.h>
#include<vtkCPProcessor.h>
#include<vtkCPPythonScriptPipeline.h>

template<typename app_type>
class ParaviewCoprocessor
{
    private:
        vtkSmartPointer<vtkCPDataDescription>   m_data_descriptor;
        vtkSmartPointer<vtkCPProcessor>         m_processor;
        int                                     m_timestep;
        int                                     m_nsteps;

    public:
        ParaviewCoprocessor()
                :       m_data_descriptor(vtkSmartPointer<vtkCPDataDescription>::New()),
                m_processor(vtkSmartPointer<vtkCPProcessor>::New()),
                m_timestep(0)
        {
            m_processor->Initialize();
            m_data_descriptor->AddInput(app().name().c_str());
        }

        app_type &app()
        {
            return *static_cast<app_type*>(this);            
        }

        ~ParaviewCoprocessor() { m_processor->Finalize(); }

        inline int initCoprocessor(int ac, char **av)
        {
            if (ac < 4)
            {
                std::cout << "Usage: " << av[0] << " <data path> <python coprocessing script> <number of time steps>" << std::endl;
                return 1;
            }
            std::string cpPythonFile = av[2];
            m_nsteps = atoi(av[3]);
            // read the coprocessing python file
            vtkSmartPointer<vtkCPPythonScriptPipeline> pipeline = vtkSmartPointer<vtkCPPythonScriptPipeline>::New();
            if (pipeline->Initialize(cpPythonFile.c_str()) == 0)
            {
                std::cout << "Problem reading the python script." << std::endl;
                return 0;
            }
            m_processor->AddPipeline(pipeline);
            return 1;
        }

        inline int run()
        {
            for (int i = 0 ; i < m_nsteps; ++i)
            {
                m_data_descriptor->SetTimeData(app().surface().time(), ++m_timestep);
                app().run();
                if (m_processor->RequestDataDescription(m_data_descriptor))
                {
                    m_data_descriptor->GetInputDescriptionByName(app().name().c_str())->SetGrid(app().vtk_storage().grid());
                    if(!m_processor->CoProcess(m_data_descriptor))
                    {
                        std::cout << "Coprocessing failed." << std::endl;
                        return 1;
                    }
                }
            }
        }

};

template<typename app_type>
int runCoprocessor(ParaviewCoprocessor<app_type> &app)
{
    return app.run();
}

#endif
