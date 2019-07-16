## Reproducible code for the paper entitled:

**Chronic fatigue syndrome and its relation with absenteeism: Elastic-net and stepwise applied to biochemical and anthropometric clinical measurements** 


Authored by:
  - [Anderson Cristiano Neisse](https://aneisse.com/)
  - [Fernando Luiz Pereira de Oliveira](http://buscatextual.cnpq.br/buscatextual/visualizacv.do?id=K4237386E4&tokenCaptchar=03AOLTBLR0gubkFsV622kUzljBk1bjblnqE22hwog9XHIZU7iuotibVz_V2TfvNRCVNQ7QBXFxl2Ob354TcvhyWUyO7RZK0-2v512U2d4LX4l0pu_pR72ieCVNhWl5PzmG5ZRNL7xz_CuonpLefzZImGERB3EyxO-nYA7zG18xcpV-iRNOOGrvKho6cFMt_ksY2lmIqyTSp3JK_NMOknjpHhgK0HbvouHiF7SfR992BBEmn6dJfhhnfsL_sw0e-28aWR2b58ZEgmC4tBrZl-XlLGMDl8jOBSNZ9-cnMcsxtvQMJAmE0lpAenvl4IsqNAH0ye4oaKkGLwhhMs0nkuvLWJ3WUbCcgypHKg)
  - [Anderson Castro Soares de Oliveira](http://buscatextual.cnpq.br/buscatextual/visualizacv.do?id=K4120983T5&tokenCaptchar=03AOLTBLS6WEn1ClxYoN8JYMmwjihJNcONdAcJBqoDHJHfxC1egMSQQqSG4C_PNlvZSU1X2NEGwKh_iX4Kkx7qNloaFuLruDYZRRS0G9iCu8b3PCT_PZvsVUGaPi7pA2AyOAnExQdreO6RAj2OKm4Lh24Egw5xUQXHMMDXLiFW_b_vLN1Ga0G1vLC7WVx00aZp9Li6cVAGZoCXdLCAjrbjiSxx9pzoTSmZ-v44lsWARwT7s1xih22-4O7yHhQj8Fb92UUqfK_t7LRqmbU6DBAUnFFV3R13Yti7Hmn7k8bz8ENTclCivJaS4oqQP_M9f5378pwKuHwxydu0u6ylUOES4ifw_Pp4JNd-3A)
  - [Frederico Rodrigues Borges da Cruz](http://buscatextual.cnpq.br/buscatextual/visualizacv.do?id=K4785689H7&tokenCaptchar=03AOLTBLTBDCE_1DYb5kchnyIYNIC7o1EYOW-jDy1BZb1jrViq5J6A2tX49GRdpHp8PY-KLJCuR7EeX3Rab6Kme9KScQ7eME_X2rwkf3OUiiGLC4tZIyXxe7gxlAw8kJVm3KHu0SMN1640T2k4CZ-mJJypQBcwmoz24ZTK7nkloE_QRMgimaMB4MhbuX7MXCcyJ7BFpICDVeus3FHUiRmCJ83AV94m5CB3NN612qVhFx9q-F_pyViovdOdF_MxxZs_ELPiSyQA-4IOjagecrSp3qZIYsts-fbFa1i5moVMqPBIMSd8EyFhEIqHD_FB4yKYPSY1mkwzEljq4Ln8cqqj84ukmozj7O3vFg)
  - [Raimundo Marques do Nascimento Neto](http://buscatextual.cnpq.br/buscatextual/visualizacv.do?id=K4705077T6&tokenCaptchar=03AOLTBLRGVQJjgQJoHiIJF_uPgO4PGOxJhhJiio5V4NRJ7wsYD606OxsTl9RUZgsGuOttOc8uu2_YgarxCXr0oVAcPTVattgw6_lvzL1iDp6ghzdOxQesvhkkX5Q__dhGV_VcUy7aHoI1B1HmSLvYBn2KGu0NmLitiDpp7W_uDEkViHThLtElI8PLL4p_Gy2ZqQeDkTN0qFaKwDLW_bo0H45xM_l9cBDcXTbSS9jFXKPGnPhRK1qFzBVaC1RTgMplIWiPcm_1fFGBtCFWbnGinIiA666sClV0FEsP3xv3J-li9BpQIOgnQFCSuBlPuNyWe5z2oUAszY2_Vwa4793QumXgtVHBgrHUSw)
  
Submitted to the *Statistical Methods in Medical Research*.

### Data availability

By request of the data provider, the data for this paper could not be shared with the code. However, some details on the dataset as well as the structure and type of variables can be found on the paper at the *Statistical Methods in Medical Research*.

### How to reproduce the results

The methodology used in the study is robust to random number generation, therefore there is no seed in the code. As a result, the final results might be slightly numerically different from the presented in the paper but there should be no difference in the conclusions.

The script `Analysis_reproduction.R` should be executed first, which sources `auxiliary_functions.R` on the run. The script `figure-table_reproduction.R` contains code to reproduce the table and figures presented in this paper.

The `roc_example.R` file is an independent file that produces the ROC curve example presented in the paper.

