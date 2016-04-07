Tina_discal_R.txt:  the R matrix from:
        
        Prob = UFget ('Pajek/Tina_DisCal') ;
        A = Prob.A ;
        p = symrcm (A'*A) ;
        A = A (:,p) ;
        R = qr (A) ;

Tina_discal_R2.txt:  the R matrix from:
        
        Prob = UFget ('Pajek/Tina_DisCal') ;
        A = Prob.A ;
        R = qr (A) ;
