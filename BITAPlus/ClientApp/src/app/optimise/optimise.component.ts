import { Component, Inject, ElementRef } from '@angular/core';
import { HttpClient, HttpHeaders } from '@angular/common/http';
import * as d3 from 'd3';
@Component({
  selector: 'app-optimise',
  templateUrl: './optimise.component.html',
  styleUrls: ['./optimise.component.css']
})
export class OptimiseComponent {
  width = 1000;//width and height for weight graph
  height = 300;
  format=d3.format('0.6f')
  opt: Optimise = {}as Optimise;
  constructor(private http: HttpClient, @Inject('BASE_URL') private baseUrl: string, public element: ElementRef) {
    http.get<Optimise>(baseUrl + 'optimise/ETL').subscribe(result => {
      this.opt=result;
      console.log(this.opt, result);
    }, error => console.error(error));
  }
  over(e: MouseEvent, inout = false) {
    d3.select(e.target as HTMLInputElement & EventTarget).classed('over', inout);
  }
  sendStep() {
    let back = +(d3.select(this.element.nativeElement).select('input.step.g').node() as HTMLInputElement & Event).value;
    console.log(back, this.opt.etLmax);
    this.opt.gamma = back;
    this.opt.tail = 0.05;
    back = +(d3.select(this.element.nativeElement).select('input.step.EU').node() as HTMLInputElement & Event).value;
    this.opt.etLmax = back;
    back = +(d3.select(this.element.nativeElement).select('input.step.EL').node() as HTMLInputElement & Event).value;
    this.opt.etLmin = back;
    this.opt.etLopt = true;
    this.sendData('optimise/ETL', this.opt)
      .subscribe(ddd => {
        console.log(ddd);
        this.opt = ddd
        console.log(this.opt.w);
      }, error => console.error(error));
  }
  sendData(key = 'optimise/ETL', sendObject: Optimise) {
    const options = {
      headers: new HttpHeaders()
        .set('Content-Type', 'application/json')
    };
    return this.http.post<Optimise>(`${this.baseUrl}${key}`, sendObject, options);
  }
}

interface Optimise {
  digit: number ,
  tdigit: number ,
  tlen: number, n: number, m: number,
  L: Array<number>,
  U: Array<number>,
  A: Array<number>,
  data: Array<number>,
  gamma: number,
  kappa: number,
  gstrength: number,
  tail: number,
  etLopt: boolean,
  etLmin: number,
  etLmax: number,
  names: Array<string>,
  w: Array<number>,
  initial: Array<number>,
  Q: Array<number>,
  etl: number,
  var: number,
  vaRindex: number,
  breakdown: Array<number>,
  back: number,
  message:string,
  mctr:Array<number>,
  risk:number,
  alpha:Array<number>,
  expreturn:number,
  CVARGLprob:boolean,
  result:Result,
  delta:number
}
interface Result{
  

  gross: number,
  longvalue: number,
  shortvalue: number,
  shortoverlong: number,
  minhold: number,
  mintrade: number,
  cost: number,
  turnover: number,
  basket: number,
  trades: number,
  var: number,
  vaRindex: number,
  etl: number,
  loss: number,
  breakdown: Array<number>,
  risk: number,
  expreturn: number,
  mctr: Array<number>,
  fmctr: Array<number>,
  sPmctr: Array<number>,
  fx: Array<number>,
  facrisk: number,
  specrisk: number
  
}