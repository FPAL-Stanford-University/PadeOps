        real(rkind), parameter :: a = (14._rkind/9._rkind)/two , b = (1._rkind/9._rkind)/four
        real(rkind) :: a06, b06
        real(rkind) :: a_np_3, b_np_3   
        real(rkind) :: a_np_2
        real(rkind) :: a_np_1, b_np_1, c_np_1, d_np_1

        a06 = a * this%onebydx; b06 = b * this%onebydx

        a_np_3 = w3*q_pp * this%onebydx  
        b_np_3 = w3*r_pp * this%onebydx 

        a_np_2 = w2*q_p * this%onebydx
        
        a_np_1 = w1*( -p * this%onebydx)
        b_np_1 = w1*(  q * this%onebydx)
        c_np_1 = w1*(  r * this%onebydx)
        d_np_1 = w1*(  s * this%onebydx)

        if (this%isBotSided) then
            rhs(:,:,1         ) =   a_np_1* fC(:,:,1         ) +  b_np_1*fC(:,:,2         )   &
                                +   c_np_1* fC(:,:,3         ) +  d_np_1*fC(:,:,4         ) 
            
            rhs(:,:,2         ) =   a_np_2*(fC(:,:,3         ) -         fC(:,:,1         ))
            
            rhs(:,:,3         ) =   a_np_3*(fC(:,:,4         ) -         fC(:,:,2         )) &
                                +   b_np_3*(fC(:,:,5         ) -         fC(:,:,1         )) 
        else
            if (this%isBotEven) then
                rhs(:,:,1) = b06*(fC(:,:,3) - fC(:,:,2)) + a06*(fC(:,:,2) - fC(:,:,1))
                rhs(:,:,2) = b06*(fC(:,:,4) - fC(:,:,1)) + a06*(fC(:,:,3) - fC(:,:,1))
                rhs(:,:,3) = b06*(fC(:,:,5) - fC(:,:,1)) + a06*(fC(:,:,4) - fC(:,:,2))
            else
                rhs(:,:,1) = b06*(fC(:,:,3) + fC(:,:,2)) + a06*(fC(:,:,2) + fC(:,:,1))
                rhs(:,:,2) = b06*(fC(:,:,4) + fC(:,:,1)) + a06*(fC(:,:,3) - fC(:,:,1))
                rhs(:,:,3) = b06*(fC(:,:,5) - fC(:,:,1)) + a06*(fC(:,:,4) - fC(:,:,2))
            end if 
        end if 

        rhs(:,:,4:this%n-3) = b06*(fC(:,:,6:this%n-1) - fC(:,:,2:this%n-5)) &
                            + a06*(fC(:,:,5:this%n-2) - fC(:,:,3:this%n-4))  

        if (this%isTopSided) then
            rhs(:,:,this%n-2  ) =   a_np_3*(fC(:,:,this%n-1  ) -         fC(:,:,this%n-3  )) &
                                +   b_np_3*(fC(:,:,this%n    ) -         fC(:,:,this%n-4  )) 
            
            rhs(:,:,this%n-1  ) =   a_np_2*(fC(:,:,this%n    ) -         fC(:,:,this%n-2  ))

            rhs(:,:,this%n    ) =  -a_np_1* fC(:,:,this%n    ) -  b_np_1*fC(:,:,this%n-1  )   &
                                -   c_np_1* fC(:,:,this%n-2  ) -  d_np_1*fC(:,:,this%n-3  )
        else
            if (this%isTopEven) then
                rhs(:,:,this%n-2) = b06*(fC(:,:,this%n  ) - fC(:,:,this%n-4)) &    
                                  + a06*(fC(:,:,this%n-1) - fC(:,:,this%n-3))    
            
                rhs(:,:,this%n-1) = b06*(fC(:,:,this%n  ) - fC(:,:,this%n-3)) &
                                  + a06*(fC(:,:,this%n  ) - fC(:,:,this%n-2))     

                rhs(:,:,this%n  ) = b06*(fC(:,:,this%n-1) - fC(:,:,this%n-2)) &
                                  + a06*(fC(:,:,this%n) - fC(:,:,this%n-1))

            else
                rhs(:,:,this%n-2) = b06*(fC(:,:,this%n  ) - fC(:,:,this%n-4)) &    
                                  + a06*(fC(:,:,this%n-1) - fC(:,:,this%n-3))    

                rhs(:,:,this%n-1) =-b06*(fC(:,:,this%n  ) + fC(:,:,this%n-3)) &
                                  + a06*(fC(:,:,this%n  ) - fC(:,:,this%n-2))     

                rhs(:,:,this%n  ) =-b06*(fC(:,:,this%n-1) + fC(:,:,this%n-2)) &
                                  - a06*(fC(:,:,this%n  ) + fC(:,:,this%n-1))
            end if 
        end if 
