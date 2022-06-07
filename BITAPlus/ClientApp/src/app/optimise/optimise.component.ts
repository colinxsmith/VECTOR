import { Component, Inject, ElementRef } from '@angular/core';
import { HttpClient, HttpHeaders } from '@angular/common/http';
import * as d3 from 'd3';
@Component({
  selector: 'app-optimise',
  templateUrl: './optimise.component.html',
  styleUrls: ['./optimise.component.css']
})
export class OptimiseComponent {
  width = 500;
  height = 400;
  opt: Array<Optimise> = [];
  constructor(private http: HttpClient, @Inject('BASE_URL') private baseUrl: string, private element: ElementRef) {
    http.get<Optimise[]>(baseUrl + 'optimise/ETL').subscribe(result => {
      this.opt.push(result[0]);
      console.log(this.opt, result);
    }, error => console.error(error));
  }
  over(e: MouseEvent, inout = false) {
    d3.select(e.target as HTMLInputElement & EventTarget).classed('over', inout);
  }
  sendStep() {
    let back = +(d3.select(this.element.nativeElement).select('input.step.g').node() as HTMLInputElement & Event).value;
    console.log(back, this.opt[0].etLmax);
    this.opt[0].gamma = back;
    back = +(d3.select(this.element.nativeElement).select('input.step.E').node() as HTMLInputElement & Event).value;
    this.opt[0].tail = 0.05;
    this.opt[0].etLmax = back;
    if (this.opt[0].etLmax != undefined) {
      this.opt[0].etLopt = true;
      this.opt[0].etLmin = back;
    } else {
      this.opt[0].etLopt = false;
    }
    this.sendData('optimise/ETL', this.opt[0])
      .subscribe(ddd => {
        console.log(ddd);
        this.opt = []
        this.opt.push(ddd[0]);
        console.log(this.opt[0].w);
      }, error => console.error(error));
  }
  sendData(key = 'optimise/ETL', sendObject: Optimise) {
    const options = {
      headers: new HttpHeaders()
        .set('Content-Type', 'application/json')
    };
    return this.http.post<Array<Optimise>>(`${this.baseUrl}${key}`, sendObject, options);
  }
}

interface Optimise {
  digit: null,
  tdigit: null,
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
  back: number
}