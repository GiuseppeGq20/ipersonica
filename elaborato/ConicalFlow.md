# Conical Flow

## Equazione di Taylor Maccoll

Studiando il problema del campo di moto inviscido attorno ad un cono a zero angolo d' attacco posto in una corrente a numero di Mach asintotico $Ma_{\infty}$  con onda d'urto attaccata al cono si giunge ad un problema differenziale ordinario del secondo ordine non lineare che è espresso dal  sistema di equazioni di Taylor Maccoll, in cui $v_r$ e $v_\omega$ sono le componenti della velocitá adimensionalizzate rispetto a $V_{lim} = \sqrt{2H}$ e $C= \frac {\gamma-1} 2 (1 - v_r^2 - v_\omega^2)$ .
$$
\begin{align}
\frac{dv_r}{d\omega}&= v_w\\
\frac{dv_\omega}{d\omega}&=\frac 1 {C - 2v_\omega²}[2v_w v_r \frac{dv_r}{d\omega} -C(2v_r + v_\omega\cot(\omega))]\\

\end{align}
$$
Questa equazioni possono essere integrate numericamente per ottenere le componenti di velocità $v_r (\omega)$ e $v_\omega (\omega)$ ,  nel codice la risoluzione di tali equazioni è implementata nella funzione `SolveTaylorMaccol() ` . La particolarità di questo sistema di equazioni sono le condizioni al contorno. Infatti per l'integrazione del sistema bisogna conoscere la velocità a valle dell'urto . Ma l' angolo d'urto è esso stesso incognito.

### Strategia risolutiva

Per risolvere il problema si puó adottare il seguente schema risolutivo:

1) si integra il sistema di Taylor Maccoll con un opportuno solutore di ODE supponendo noto l'angolo d'urto $\beta$ e l'integrazione si ferma al soddisfacimento della condizione $v_w=0$ , ovvero la condizione di tangenza della corrente sul corpo. Ovviamente l'angolo $\omega_{fin}$ per cui si ha  $v_w=0$ non coincide con l'angolo $\delta_c$ del cono.
   Nel codice tale routine é implementata nella funzione `SolveTaylorMaccol()`

2) Si può a questo punto utilizzare un algoritmo di root-finding sulla funzione errore:
   $$
   e(\beta)= \omega_{fin}(\beta)-\delta_c
   $$
   l'angolo $\beta$ cosí ottenuto é l'angolo dell'onda d'urto conica.
   Nel codice questo é implementato con un metodo delle secanti nella funzione `betaCone()`, al variare dell'inizializzazione dell' algoritmo si possono ottenere la soluzione forte o debole dell'urto, tuttavia il metodo non garantisce convergenza per ogni valore di tentativo per l'inizializzazione del metodo delle secanti.

   

## massimo angolo di semiapertura del cono a fissato $Ma_\infty$

Il massimo angolo di apertura del cono $\delta_c$ corrisponde a quell' angolo per cui il numero numero di mach a valle dell' onda d' urto è unitario. 
Si può impostare la seguente strategia per risolvere questo problema:

- definire la funzione $f(Ma_\infty , \beta) = M_2 -1$
- trovare l' angolo $\beta$ per cui $f=0$ 
- Integrare le equazioni di Taylor-Maccoll, noti $\beta$ e $Ma_\infty$ , da cui si può ottenere $\delta_{c,max}$

## calcolo coefficienti aerodinamici

Attraverso i metodo del cono locale e di High é possibile calcolare il $c_p$ sulla superficie del cono. Per la genesi di questi metodi, il $c_p$ e funzione solo dell'angolo $\phi$ , da questa osservazione si ha che la forza aerodinamica per unità di pressione dinamica é data dalla relazione $\ref{eq:Cp}$:
$$
\frac F {q_\infty}= \int_0^R \int_0^{2\pi}  c_p(\phi) \frac r {sin(\delta_c)}drd\phi = \frac {R^2} {2 sin(\delta_c)} \int_0^{2\pi} c_p(\phi) d\phi 		\tag{1} \label{eq:Cp}
$$

### calcolo $C_{Fx}$ e ${C_{Fy}}$

L' asse x coincide con l'asse del cono e dunque la forza per unitá di pressione dinamica lungo l'asse x,  $F_x$,  é data dalla relazione $\ref{Fx}$ .
Per ottenere la forza per unitá di pressione dinamica lungo l'asse y , $F_y$, bisogna proiettare la forza aerodinamica radiale $F_R = Fcos(\delta_c)$ lungo l'asse y ottenendo la relazione $\ref{Fy}$.
$$
\begin{equation} \label{Fx} \tag{2}
C_{F_x}= \frac {F_{x}}{ \frac {\pi R^2}{\sin(\delta_c)}} = \frac {\sin(\delta_c)}{2\pi} \int_0^{2\pi} c_pd\phi
\end{equation}
$$

$$
\begin{equation} \label{Fy} \tag{3}
C_{F_y}= \frac {F_{y}}{ \frac {\pi R^2}{\sin(\delta_c)}} = \frac {\cos(\delta_c)}{2\pi} \int_0^{2\pi} c_p\cos(\phi)d\phi
\end{equation}
$$


Questi due integrali posso essere calcolati numericamente una volta nota la distribuzione di $c_p$. 

### calcolo $C_l$ e $C_d$

I coefficienti $C_l$ e $C_d$ sono calcolati attraverso la seguente trasformazione di rotazione di un angolo  pari all'angolo di attacco $\alpha$ :
$$
\begin{align}
C_d&= C_{F_x}\cos(\alpha) + C_{F_y}\sin(\alpha) \\
C_l&= -C_{F_x}\sin(\alpha) + C_{F_y}\cos(\alpha)\\
\end{align}
$$
  Il calcolo dei coefficienti aerodinamici é implementato nella funzione `calcClCdCone()`. Gli integrali che compaiono nell eqs. $\ref{Fx}$ e $\ref{Fy}$ sono calcolati attraverso una regola trapezoidale. 



