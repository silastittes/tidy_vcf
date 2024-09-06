from setuptools import setup, find_packages 
  
long_description = 'Given a vcf file, produces a tidy versions of sites and genotypes data, in \
                    efforts to make it easier to calculate summary statitics and visualizations.' 
  
setup( 
        name ='tidy-vcf', 
        version ='0.2.1.0', 
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
                'tidy-vcf=tidy_vcf.tidy_vcf:main'
            ] 
        }, 
        test_suite="tests",
        tests_require=["pytest"],
        classifiers =( 
            "Programming Language :: Python :: 3", 
            "License :: OSI Approved :: MIT License", 
            "Operating System :: OS Independent", 
        ), 
        keywords ='VCF genetics python package', 
        zip_safe = False,
        python_requires='>=3.7'
) 
