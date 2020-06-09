from setuptools import setup, find_packages 
  
with open('requirements.txt') as f: 
    requirements = f.readlines() 
  
long_description = 'Given a vcf file, produces a tidy versions of sites and genotypes data, in \
                    efforts to make it easier to calculate summary statitics and visualizations.' 
  
setup( 
        name ='tidy_vcf', 
        version ='0.1.2.3', 
        author ='Silas Tittes', 
        author_email ='silas.tittes@gmail.com', 
        url ='https://github.com/silastittes/tidy_vcf', 
        description ='Make tidy VCF data.', 
        long_description = long_description, 
        long_description_content_type ="text/markdown", 
        license ='MIT', 
        packages = find_packages(), 
        entry_points ={ 
            'console_scripts': [ 
                'tidy_vcf = tidy_vcf.tidy_vcf:main'
            ] 
        }, 
        classifiers =( 
            "Programming Language :: Python :: 3", 
            "License :: OSI Approved :: MIT License", 
            "Operating System :: OS Independent", 
        ), 
        keywords ='VCF genetics python package', 
        install_requires = requirements, 
        zip_safe = False,
        python_requires='>=3.6'
) 
